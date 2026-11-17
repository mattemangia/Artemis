using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Artemis.Particles
{
    // CPU solver of the massive particle system.
    //
    // The solver is position based: each sub-step predicts positions from the
    // velocities, pushes overlapping particles apart, clamps them to the bounds and
    // then derives the new velocities from the distance actually travelled. Every
    // particle computes its own correction from a read-only snapshot of its
    // neighbours (a Jacobi iteration), so the work can be spread over all cores
    // without two threads ever writing to the same particle, and the result does not
    // depend on thread scheduling.
    //
    // Two things keep very large scenes affordable:
    //   * frozen particles (a mesh turned into sand that has not been released yet)
    //     and sleeping particles (sand that has come to rest) are skipped by the
    //     contact solver and only act as obstacles;
    //   * neighbours are found through a hashed uniform grid built with a counting
    //     sort, whose memory use depends on the particle count and not on the size of
    //     the world.
    public partial class MassiveParticleSystem
    {
        #region Flags

        private const byte FlagAlive = 1;
        private const byte FlagCollides = 2;
        private const byte FlagGravity = 4;
        private const byte FlagFrozen = 8;
        private const byte FlagSleeping = 16;

        private const byte WakeUp = 1;
        private const byte WakeAndRelease = 2;

        #endregion

        #region Solver Fields

        // One record per colliding particle, stored in cell order. Everything the
        // neighbour loop needs sits together, so testing a candidate touches a single
        // cache line.
        private struct SortedParticle
        {
            public float X, Y, Z;
            public float Radius;
            public float Weight;   // inverse mass, or -1 frozen / -2 asleep
            public float MoveX, MoveY, MoveZ;
            public int Key;
            public int Index;
            public int Group;      // rigid group, or -1
        }

        // Per-particle solver state.
        private float[] _dpX = Array.Empty<float>();
        private float[] _dpY = Array.Empty<float>();
        private float[] _dpZ = Array.Empty<float>();
        private byte[] _rest = Array.Empty<byte>();
        private byte[] _wake = Array.Empty<byte>();
        private byte[] _hit = Array.Empty<byte>();

        // Particles sorted by grid cell: a compact copy of what the neighbour search
        // reads, so the inner loop walks contiguous memory.
        private SortedParticle[] _sorted = Array.Empty<SortedParticle>();
        private float[] _stabX = Array.Empty<float>();
        private float[] _stabY = Array.Empty<float>();
        private float[] _stabZ = Array.Empty<float>();

        private int[] _bucketStart = Array.Empty<int>();
        private int[] _bucketCursor = Array.Empty<int>();
        private int[][] _threadHistograms = { Array.Empty<int>() };
        private int _bucketMask;
        private bool _gridDense;
        private int _gridX0, _gridY0, _gridZ0, _gridNx, _gridNy;
        private int _sortedCount;
        private float _gridInvCell;
        private float _maxRadius;
        private int _activeCount;
        private int _awakeColliders;

        #endregion

        #region Solver Properties

        /// <summary>
        /// Whether particles that have come to rest are put to sleep. Sleeping
        /// particles are skipped by the solver until something disturbs them.
        /// </summary>
        public bool AllowSleeping { get; set; } = true;

        /// <summary>
        /// Speed (m/s) below which a particle counts as resting.
        /// </summary>
        public float SleepSpeed { get; set; } = 0.05f;

        /// <summary>
        /// Speed (m/s) a particle needs in order to wake sleeping neighbours. It is
        /// higher than <see cref="SleepSpeed"/> so that the faint jitter of grains
        /// settling next to a sleeping pile does not keep the whole pile awake.
        /// </summary>
        public float WakeSpeed { get; set; } = 0.1f;

        /// <summary>
        /// Number of consecutive resting sub-steps after which a particle falls asleep.
        /// </summary>
        public int SleepSubSteps { get; set; } = 120;

        /// <summary>
        /// Impact speed (m/s) at which a moving particle knocks a frozen particle
        /// loose. The default (infinity) means frozen particles are only released
        /// explicitly, for example with <see cref="ReleaseAll"/>.
        /// </summary>
        public float FrozenBreakSpeed { get; set; } = float.PositiveInfinity;

        /// <summary>
        /// Fraction of each overlap removed per sub-step. Lower values are softer but
        /// calmer; values above 1 over-correct and can jitter.
        /// </summary>
        public float ContactRelaxation { get; set; } = 1.0f;

        /// <summary>
        /// How strongly a contact between two grains moves the upper one rather than
        /// the lower one (0 = split evenly, 1 = the upper grain does all the moving
        /// when it sits directly on top). Without this bias the weight of a tall pile
        /// compresses its bottom layers.
        /// </summary>
        public float StackingBias { get; set; } = 0.9f;

        /// <summary>
        /// Number of particles the solver actually moved in the last sub-step, i.e.
        /// alive particles that are neither frozen nor asleep.
        /// </summary>
        public int ActiveCount => _activeCount;

        /// <summary>
        /// Number of frozen particles.
        /// </summary>
        public int FrozenCount => CountFlag(FlagFrozen);

        /// <summary>
        /// Number of sleeping particles.
        /// </summary>
        public int SleepingCount => CountFlag(FlagSleeping);

        #endregion

        #region Allocation

        private void AllocateSolverArrays(int capacity)
        {
            _dpX = new float[capacity];
            _dpY = new float[capacity];
            _dpZ = new float[capacity];
            _rest = new byte[capacity];
            _wake = new byte[capacity];
            _hit = new byte[capacity];
            AllocateGroupArrays(capacity);
            AllocateSoftArrays(capacity);

            _sorted = new SortedParticle[capacity];
            _stabX = new float[capacity];
            _stabY = new float[capacity];
            _stabZ = new float[capacity];
        }

        private int CountFlag(byte flag)
        {
            SyncHost();
            int total = 0;
            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & (FlagAlive | flag)) == (FlagAlive | flag))
                    total++;
            }
            return total;
        }

        #endregion

        #region Frozen Particles

        /// <summary>
        /// Freezes or releases a single particle. A frozen particle keeps its place
        /// and acts as an immovable obstacle until it is released.
        /// </summary>
        public void SetFrozen(int index, bool frozen)
        {
            TouchHost();
            if (index < 0 || index >= _count)
                throw new ArgumentOutOfRangeException(nameof(index));

            if (frozen)
            {
                _flags[index] |= FlagFrozen;
                _velX[index] = _velY[index] = _velZ[index] = 0;
            }
            else
            {
                _flags[index] &= unchecked((byte)~(FlagFrozen | FlagSleeping));
                _rest[index] = 0;
            }
        }

        /// <summary>
        /// Freezes a range of particles, for example the ones just spawned from a mesh.
        /// </summary>
        public void Freeze(int startIndex, int count)
        {
            TouchHost();
            int end = Math.Min(startIndex + count, _count);
            for (int i = Math.Max(startIndex, 0); i < end; i++)
            {
                _flags[i] |= FlagFrozen;
                _velX[i] = _velY[i] = _velZ[i] = 0;
            }
        }

        /// <summary>
        /// Releases every frozen particle: the shape collapses into loose sand.
        /// </summary>
        /// <returns>The number of particles released.</returns>
        public int ReleaseAll() => ReleaseWhere(static (_, _, _) => true);

        /// <summary>
        /// Releases the frozen particles inside a sphere.
        /// </summary>
        /// <returns>The number of particles released.</returns>
        public int ReleaseInRadius(Vector3 center, float radius)
        {
            float radiusSq = radius * radius;
            return ReleaseWhere((x, y, z) =>
            {
                float dx = x - center.X, dy = y - center.Y, dz = z - center.Z;
                return dx * dx + dy * dy + dz * dz <= radiusSq;
            });
        }

        /// <summary>
        /// Releases the frozen particles on the positive side of a plane
        /// (<c>dot(normal, p) &gt; offset</c>). Sweeping the offset over time makes a
        /// shape crumble progressively, for example from the top down.
        /// </summary>
        /// <returns>The number of particles released.</returns>
        public int ReleaseBeyondPlane(Vector3 normal, float offset)
            => ReleaseWhere((x, y, z) => normal.X * x + normal.Y * y + normal.Z * z > offset);

        /// <summary>
        /// Releases the frozen particles whose position satisfies a predicate.
        /// </summary>
        /// <returns>The number of particles released.</returns>
        public int ReleaseWhere(Func<float, float, float, bool> predicate)
        {
            TouchHost();
            if (predicate == null) throw new ArgumentNullException(nameof(predicate));

            int released = 0;
            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & (FlagAlive | FlagFrozen)) != (FlagAlive | FlagFrozen))
                    continue;
                if (!predicate(_posX[i], _posY[i], _posZ[i]))
                    continue;

                _flags[i] &= unchecked((byte)~(FlagFrozen | FlagSleeping));
                _rest[i] = 0;
                released++;
            }
            return released;
        }

        /// <summary>
        /// Wakes every sleeping particle (for example after changing gravity or bounds).
        /// </summary>
        public void WakeAll()
        {
            TouchHost();
            for (int i = 0; i < _count; i++)
            {
                _flags[i] &= unchecked((byte)~FlagSleeping);
                _rest[i] = 0;
            }
        }

        /// <summary>
        /// Whether a particle is frozen.
        /// </summary>
        public bool IsFrozen(int index) { SyncHost(); return (_flags[index] & FlagFrozen) != 0; }

        /// <summary>
        /// Whether a particle is asleep.
        /// </summary>
        public bool IsSleeping(int index) { SyncHost(); return (_flags[index] & FlagSleeping) != 0; }

        #endregion

        #region Sub-step


        private void UpdateCpu(float dt)
        {
            if (dt <= 0)
                return;

            Predict(dt);

            // Shape matching first, contacts after: the contacts get the last word, so
            // a soft body cannot pull its particles back into whatever it touches.
            MatchSoftBodies(dt);

            // With every particle frozen or asleep there is nothing to solve: an idle
            // scene costs one pass over the flags.
            if (EnableCollisions && _maxRadius > 0 && _awakeColliders > 0)
            {
                BuildGrid();
                SolveContacts(dt);
            }

            MatchRigidGroups(dt);
            Finish(dt);
            DampSoftBodies();
        }

        private void ForRanges(int count, Action<int, int> body)
        {
            if (count <= 0)
                return;

            const int minChunk = 2048;
            if (count <= minChunk || _threadCount <= 1)
            {
                body(0, count);
                return;
            }

            // Several chunks per thread so an uneven load (a dense pile next to
            // sleeping sand) is balanced by work stealing.
            int chunk = Math.Max(minChunk, count / (_threadCount * 8));
            Parallel.ForEach(Partitioner.Create(0, count, chunk), _parallelOptions,
                range => body(range.Item1, range.Item2));
        }

        /// <summary>
        /// Applies pending wake-ups, integrates forces and predicts the new positions.
        /// </summary>
        private void Predict(float dt)
        {
            float gx = Gravity.X * dt, gy = Gravity.Y * dt, gz = Gravity.Z * dt;
            float damping = Damping;
            int awakeColliders = 0;

            ForRanges(_count, (start, end) =>
            {
                int localAwake = 0;

                for (int i = start; i < end; i++)
                {
                    byte flags = _flags[i];
                    _hit[i] = 0;

                    byte wake = _wake[i];
                    if (wake != 0)
                    {
                        _wake[i] = 0;
                        flags &= unchecked((byte)~FlagSleeping);
                        if (wake == WakeAndRelease)
                            flags &= unchecked((byte)~FlagFrozen);
                        _flags[i] = flags;
                        _rest[i] = 0;
                    }

                    _prevX[i] = _posX[i];
                    _prevY[i] = _posY[i];
                    _prevZ[i] = _posZ[i];

                    if ((flags & FlagAlive) == 0 || (flags & (FlagFrozen | FlagSleeping)) != 0 || _invMass[i] <= 0)
                        continue;

                    if ((flags & FlagCollides) != 0)
                        localAwake++;

                    float vx = _velX[i], vy = _velY[i], vz = _velZ[i];
                    if ((flags & FlagGravity) != 0)
                    {
                        vx += gx;
                        vy += gy;
                        vz += gz;
                    }

                    vx *= damping;
                    vy *= damping;
                    vz *= damping;

                    _velX[i] = vx;
                    _velY[i] = vy;
                    _velZ[i] = vz;
                    _posX[i] += vx * dt;
                    _posY[i] += vy * dt;
                    _posZ[i] += vz * dt;
                }

                System.Threading.Interlocked.Add(ref awakeColliders, localAwake);
            });

            _awakeColliders = awakeColliders;
        }

        #endregion

        #region Hashed Grid

        // Rows of cells (fixed y and z) are hashed; the cells of a row occupy
        // consecutive buckets. The three cells x-1, x, x+1 of a row are therefore one
        // contiguous range of the sorted arrays, which cuts the neighbour search from
        // 27 table lookups to 9.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int HashRow(int cy, int cz)
            => unchecked((cy * 19349663) ^ (cz * 83492791));

        // Low bits of the cell coordinates. Unrelated cells can land in the same
        // bucket; the key tells them apart from the cell that was asked for (the 27
        // cells around a particle always have distinct keys).
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int CellKey(int cx, int cy, int cz)
            => (cx & 1023) | ((cy & 1023) << 10) | ((cz & 1023) << 20);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int CellCoord(float value, float invCell)
            => (int)MathF.Floor(value * invCell);

        /// <summary>
        /// Sorts the colliding particles by grid cell (counting sort over a hash table)
        /// and copies what the contact solver needs into cell order.
        /// </summary>
        private void BuildGrid()
        {
            int count = _count;

            // Cells are one particle diameter wide, so a particle only ever touches
            // particles in its own and the 26 surrounding cells.
            float cell = Math.Max(2f * _maxRadius, 1e-4f);
            float invCell = 1f / cell;
            _gridInvCell = invCell;

            // Extent of the colliding particles, in cells.
            int minCx = int.MaxValue, minCy = int.MaxValue, minCz = int.MaxValue;
            int maxCx = int.MinValue, maxCy = int.MinValue, maxCz = int.MinValue;
            object gate = new();
            ForRanges(count, (start, end) =>
            {
                float lx = float.MaxValue, ly = float.MaxValue, lz = float.MaxValue;
                float hx = float.MinValue, hy = float.MinValue, hz = float.MinValue;
                for (int i = start; i < end; i++)
                {
                    if ((_flags[i] & (FlagAlive | FlagCollides)) != (FlagAlive | FlagCollides))
                        continue;

                    float x = _posX[i], y = _posY[i], z = _posZ[i];
                    if (x < lx) lx = x; if (x > hx) hx = x;
                    if (y < ly) ly = y; if (y > hy) hy = y;
                    if (z < lz) lz = z; if (z > hz) hz = z;
                }

                if (lx > hx)
                    return;

                lock (gate)
                {
                    minCx = Math.Min(minCx, CellCoord(lx, invCell)); maxCx = Math.Max(maxCx, CellCoord(hx, invCell));
                    minCy = Math.Min(minCy, CellCoord(ly, invCell)); maxCy = Math.Max(maxCy, CellCoord(hy, invCell));
                    minCz = Math.Min(minCz, CellCoord(lz, invCell)); maxCz = Math.Max(maxCz, CellCoord(hz, invCell));
                }
            });

            // Preferred layout: a plain linear grid over the occupied extent, padded
            // by one cell so that neighbour lookups never leave it. Records then end
            // up in spatial order, and while the solver sweeps through them the nine
            // rows it reads advance through memory in step - which the CPU prefetches
            // well. Only when the particles are spread so thinly that such a grid
            // would be mostly empty is a hash table used instead.
            long nx = (long)maxCx - minCx + 3, ny = (long)maxCy - minCy + 3, nz = (long)maxCz - minCz + 3;
            long denseCells = minCx > maxCx ? 0 : nx * ny * nz;
            bool dense = denseCells > 0 && denseCells <= Math.Max(8L * count, 1 << 16);

            int tableSize;
            if (dense)
            {
                tableSize = (int)denseCells;
                _gridX0 = minCx - 1;
                _gridY0 = minCy - 1;
                _gridZ0 = minCz - 1;
                _gridNx = (int)nx;
                _gridNy = (int)ny;
            }
            else
            {
                tableSize = 1024;
                while (tableSize < count * 2 && tableSize < (1 << 28))
                    tableSize <<= 1;
            }
            _gridDense = dense;

            // Three extra entries let a row be read as one run of three buckets.
            if (_bucketStart.Length < tableSize + 4)
            {
                _bucketStart = new int[tableSize + 4];
                _bucketCursor = new int[tableSize + 4];
            }
            else
            {
                Array.Clear(_bucketStart, 0, tableSize + 4);
            }

            int mask = tableSize - 1;
            _bucketMask = mask;
            int x0 = _gridX0, y0 = _gridY0, z0 = _gridZ0, gnx = _gridNx, gny = _gridNy;

            // Bucket of every particle (-1 for the ones that take no part in collisions).
            ForRanges(count, (start, end) =>
            {
                for (int i = start; i < end; i++)
                {
                    if ((_flags[i] & (FlagAlive | FlagCollides)) != (FlagAlive | FlagCollides))
                    {
                        _cellIndices[i] = -1;
                        continue;
                    }

                    int cx = CellCoord(_posX[i], invCell);
                    int cy = CellCoord(_posY[i], invCell);
                    int cz = CellCoord(_posZ[i], invCell);
                    _cellIndices[i] = dense
                        ? ((cz - z0) * gny + (cy - y0)) * gnx + (cx - x0)
                        : unchecked(HashRow(cy, cz) + cx) & mask;
                }
            });

            if (count >= 100_000 && _threadCount > 1)
                SortParallel(count, tableSize);
            else
                SortSequential(count, tableSize);

            // Gather the data the neighbour loop reads.
            ForRanges(_sortedCount, (start, end) =>
            {
                var sorted = _sorted;
                for (int s = start; s < end; s++)
                {
                    ref var record = ref sorted[s];
                    int i = record.Index;
                    byte flags = _flags[i];
                    bool movable = (flags & (FlagFrozen | FlagSleeping)) == 0;

                    float x = _posX[i], y = _posY[i], z = _posZ[i];
                    record.X = x;
                    record.Y = y;
                    record.Z = z;
                    record.Radius = _radius[i];
                    record.Key = CellKey(CellCoord(x, invCell), CellCoord(y, invCell), CellCoord(z, invCell));

                    // Negative weights tag obstacles: -1 frozen, -2 asleep.
                    record.Weight = movable
                        ? _invMass[i]
                        : (flags & FlagFrozen) != 0 ? -1f : -2f;

                    record.Group = BodyTag(i);
                    record.MoveX = x - _prevX[i];
                    record.MoveY = y - _prevY[i];
                    record.MoveZ = z - _prevZ[i];
                }
            });
        }

        private void SortSequential(int count, int tableSize)
        {
            // Histogram and prefix sum.
            var bucketStart = _bucketStart;
            for (int i = 0; i < count; i++)
            {
                int bucket = _cellIndices[i];
                if (bucket >= 0)
                    bucketStart[bucket + 1]++;
            }

            for (int b = 0; b < tableSize; b++)
                bucketStart[b + 1] += bucketStart[b];

            _sortedCount = bucketStart[tableSize];
            bucketStart[tableSize + 1] = bucketStart[tableSize + 2] = bucketStart[tableSize + 3] = _sortedCount;

            // Scatter particle indices into cell order.
            var cursor = _bucketCursor;
            Array.Copy(bucketStart, cursor, tableSize);
            for (int i = 0; i < count; i++)
            {
                int bucket = _cellIndices[i];
                if (bucket >= 0)
                    _sorted[cursor[bucket]++].Index = i;
            }
        }

        /// <summary>
        /// Counting sort with one histogram per thread. Each thread owns a contiguous
        /// slice of the particles; within a bucket the slices are laid out in thread
        /// order, so the result is identical to the sequential sort.
        /// </summary>
        private void SortParallel(int count, int tableSize)
        {
            int threads = Math.Min(_threadCount, 16);
            if (_threadHistograms.Length != threads || _threadHistograms[0].Length < tableSize)
            {
                _threadHistograms = new int[threads][];
                for (int t = 0; t < threads; t++)
                    _threadHistograms[t] = new int[tableSize];
            }

            var histograms = _threadHistograms;

            Parallel.For(0, threads, _parallelOptions, t =>
            {
                var histogram = histograms[t];
                Array.Clear(histogram, 0, tableSize);

                int from = (int)((long)count * t / threads);
                int to = (int)((long)count * (t + 1) / threads);
                for (int i = from; i < to; i++)
                {
                    int bucket = _cellIndices[i];
                    if (bucket >= 0)
                        histogram[bucket]++;
                }
            });

            // Turn the per-thread counts into per-thread write cursors.
            var bucketStart = _bucketStart;
            int running = 0;
            for (int b = 0; b < tableSize; b++)
            {
                bucketStart[b] = running;
                for (int t = 0; t < threads; t++)
                {
                    int n = histograms[t][b];
                    histograms[t][b] = running;
                    running += n;
                }
            }
            bucketStart[tableSize] = running;
            bucketStart[tableSize + 1] = bucketStart[tableSize + 2] = bucketStart[tableSize + 3] = running;
            _sortedCount = running;

            Parallel.For(0, threads, _parallelOptions, t =>
            {
                var cursor = histograms[t];
                var sorted = _sorted;

                int from = (int)((long)count * t / threads);
                int to = (int)((long)count * (t + 1) / threads);
                for (int i = from; i < to; i++)
                {
                    int bucket = _cellIndices[i];
                    if (bucket >= 0)
                        sorted[cursor[bucket]++].Index = i;
                }
            });
        }

        #endregion

        #region Contacts

        /// <summary>
        /// Computes, for every moving particle, the displacement that resolves its
        /// overlaps (with Coulomb friction) and applies it.
        /// </summary>
        private void SolveContacts(float dt)
        {
            ForRanges(_sortedCount, (start, end) => SolveRange(start, end, dt));

            ForRanges(_sortedCount, (start, end) =>
            {
                for (int s = start; s < end; s++)
                {
                    if (_sorted[s].Weight <= 0)
                        continue;

                    int i = _sorted[s].Index;
                    _posX[i] += _dpX[i];
                    _posY[i] += _dpY[i];
                    _posZ[i] += _dpZ[i];
                }
            });
        }

        /// <summary>
        /// Contact pass for a range of sorted records. Kept as a method (rather than a
        /// lambda body) so that all parameters are plain locals the JIT can keep in
        /// registers instead of fields of a closure.
        /// </summary>
        [MethodImpl(MethodImplOptionsCompat.AggressiveOptimization)]
        private void SolveRange(int start, int end, float dt)
        {
            GetContactParameters(dt, out float friction, out float relaxation, out float movingThresholdSq,
                out float breakDistance, out bool canBreak, out float stackingBias, out float upX, out float upY, out float upZ);
            float invCell = _gridInvCell;
            int mask = _bucketMask;
            bool dense = _gridDense;
            int x0 = _gridX0, y0 = _gridY0, z0 = _gridZ0, gnx = _gridNx, gny = _gridNy;
            int[] bucketStart = _bucketStart;
            SortedParticle[] sorted = _sorted;
            byte[] wakeFlags = _wake;
            float[] outX = _dpX, outY = _dpY, outZ = _dpZ;
            float[] outStabX = _stabX, outStabY = _stabY, outStabZ = _stabZ;

            for (int s = start; s < end; s++)
            {
                ref readonly var self = ref sorted[s];
                float wi = self.Weight;
                int selfGroup = self.Group;
                if (wi <= 0)
                    continue; // frozen, asleep or infinitely heavy: nothing to solve

                float px = self.X, py = self.Y, pz = self.Z;
                float ri = self.Radius;
                float mx = self.MoveX, my = self.MoveY, mz = self.MoveZ;
                bool moving = mx * mx + my * my + mz * mz > movingThresholdSq;

                int cx = CellCoord(px, invCell);
                int cy = CellCoord(py, invCell);
                int cz = CellCoord(pz, invCell);

                float dpx = 0, dpy = 0, dpz = 0;
                float stx = 0, sty = 0, stz = 0;
                int contacts = 0;
                float nxx = 0, nxy = 0, nxz = 0, nyy = 0, nyz = 0, nzz = 0;
                float fx = 0, fy = 0, fz = 0;

                for (int oz = -1; oz <= 1; oz++)
                {
                    for (int oy = -1; oy <= 1; oy++)
                    {
                        int rowHash = dense ? 0 : HashRow(cy + oy, cz + oz);
                        int rowKey = (((cy + oy) & 1023) << 10) | (((cz + oz) & 1023) << 20);

                        // The cells x-1, x, x+1 of a row are consecutive buckets,
                        // i.e. one contiguous run of records - unless the run
                        // wraps around the end of the table, in which case the
                        // row is walked in up to three pieces.
                        int first, pieces;
                        if (dense)
                        {
                            first = ((cz + oz - z0) * gny + (cy + oy - y0)) * gnx + (cx - 1 - x0);
                            pieces = 1;
                        }
                        else
                        {
                            first = unchecked(rowHash + cx - 1) & mask;
                            pieces = first <= mask - 2 ? 1 : 3;
                        }

                        for (int piece = 0; piece < pieces; piece++)
                        {
                            int from, to;
                            if (pieces == 1)
                            {
                                from = bucketStart[first];
                                to = bucketStart[first + 3];
                            }
                            else
                            {
                                int bucket = (first + piece) & mask;
                                from = bucketStart[bucket];
                                to = bucketStart[bucket + 1];
                            }

                            var run = new ReadOnlySpan<SortedParticle>(sorted, from, to - from);
                            for (int u = 0; u < run.Length; u++)
                            {
                                ref readonly var other = ref run[u];
                                int t = from + u;

                                // Reject records of unrelated cells that share
                                // these buckets: same row, and x within one cell.
                                int key = other.Key;

                                // Particles of one rigid object do not push each other.
                                if (selfGroup >= 0 && other.Group == selfGroup)
                                    continue;

                                if (t == s || (!dense && (((key ^ rowKey) & ~1023) != 0 || ((key - cx + 1) & 1023) > 2)))
                                    continue;

                                float dx = px - other.X;
                                float dy = py - other.Y;
                                float dz = pz - other.Z;
                                float minDist = ri + other.Radius;
                                float distSq = dx * dx + dy * dy + dz * dz;
                                float wj = other.Weight;

                                // Something moving right next to sleeping sand wakes
                                // it (with a small margin, since resting grains
                                // barely touch), so that nothing stays asleep in
                                // mid-air when its support is carried away.
                                if (wj == -2f && moving && distSq < minDist * minDist * 1.21f)
                                    { if (wakeFlags[other.Index] < WakeUp) wakeFlags[other.Index] = WakeUp; }

                                if (distSq >= minDist * minDist)
                                    continue;

                                float nx, ny, nz, dist;
                                if (distSq > 1e-12f)
                                {
                                    dist = MathF.Sqrt(distSq);
                                    float inv = 1f / dist;
                                    nx = dx * inv; ny = dy * inv; nz = dz * inv;
                                }
                                else
                                {
                                    // Coincident centres: separate along a fixed
                                    // axis, in opposite directions for the pair.
                                    dist = 0;
                                    nx = 0; ny = self.Index > other.Index ? 1f : -1f; nz = 0;
                                }

                                // Relative movement during this sub-step.
                                float rx = mx, ry = my, rz = mz;
                                if (wj > 0)
                                {
                                    rx -= other.MoveX;
                                    ry -= other.MoveY;
                                    rz -= other.MoveZ;
                                }
                                float rn = rx * nx + ry * ny + rz * nz;

                                // A fast grain knocks a frozen one loose.
                                if (wj == -1f && canBreak && -rn > breakDistance)
                                    wakeFlags[other.Index] = WakeAndRelease;

                                float penetration = minDist - dist;
                                float wOther = wj > 0 ? wj : 0;
                                float share = wi / (wi + wOther);

                                // In a pile the weight of everything above has to be
                                // carried by the grains below. Splitting each overlap
                                // evenly lets that load squash the bottom layers, so
                                // the grain that sits higher takes most of the
                                // correction, as if the supporting grain were heavier.
                                if (wOther > 0 && stackingBias > 0)
                                {
                                    float up = nx * upX + ny * upY + nz * upZ;
                                    share = Math.Clamp(share * (1f + stackingBias * up), 0f, 1f);
                                }

                                float magnitude = penetration * share;
                                dpx += nx * magnitude;
                                dpy += ny * magnitude;
                                dpz += nz * magnitude;
                                nxx += nx * nx; nxy += nx * ny; nxz += nx * nz;
                                nyy += ny * ny; nyz += ny * nz; nzz += nz * nz;
                                contacts++;

                                // Only the overlap created by this sub-step's approach
                                // may turn into velocity. Overlap that was already
                                // there (spawned intersecting, or left over from an
                                // unconverged pile) is removed without a kick,
                                // otherwise dense packs would blow themselves apart.
                                float approach = rn < 0 ? -rn : 0;
                                float preexisting = penetration - approach;
                                if (preexisting > 0)
                                {
                                    stx += nx * preexisting * share;
                                    sty += ny * preexisting * share;
                                    stz += nz * preexisting * share;
                                }

                                if (friction > 0)
                                {
                                    float tx = rx - rn * nx, ty = ry - rn * ny, tz = rz - rn * nz;
                                    float slideSq = tx * tx + ty * ty + tz * tz;
                                    if (slideSq > 1e-16f)
                                    {
                                        // Coulomb: friction can undo at most
                                        // mu * penetration of the slide.
                                        float slide = MathF.Sqrt(slideSq);
                                        float limit = friction * penetration;
                                        float scale = slide <= limit ? 1f : limit / slide;
                                        scale *= share;

                                        fx -= tx * scale;
                                        fy -= ty * scale;
                                        fz -= tz * scale;
                                    }
                                }
                            }
                        }
                    }
                }

                int index = self.Index;
                if (contacts > 0)
                {
                    // Every neighbour asks for its own overlap to be removed along
                    // its own normal. Simply adding the requests overshoots (three
                    // supports would lift a grain twice as far as needed) and
                    // averaging them undershoots (grains sink). The step along the
                    // summed direction D that satisfies the requests best in the
                    // least-squares sense is alpha * D, with
                    //   alpha = |D|^2 / sum((n . D)^2).
                    // It is exact for a single contact and for symmetric supports.
                    float lengthSq = dpx * dpx + dpy * dpy + dpz * dpz;
                    float projected =
                        dpx * (nxx * dpx + nxy * dpy + nxz * dpz) +
                        dpy * (nxy * dpx + nyy * dpy + nyz * dpz) +
                        dpz * (nxz * dpx + nyz * dpy + nzz * dpz);
                    float alpha = projected > 1e-20f ? Math.Min(lengthSq / projected, 1.5f) : 0f;
                    float scale = relaxation * alpha;

                    // Friction: the contacts share the job of holding the grain, so
                    // together they cannot undo more than the grain actually moved.
                    float frictionScale = Math.Min(alpha, 1f);
                    fx *= frictionScale; fy *= frictionScale; fz *= frictionScale;
                    float frictionSq = fx * fx + fy * fy + fz * fz;
                    float movedSq = mx * mx + my * my + mz * mz;
                    if (frictionSq > movedSq && frictionSq > 0)
                    {
                        float clamp = MathF.Sqrt(movedSq / frictionSq);
                        fx *= clamp; fy *= clamp; fz *= clamp;
                    }

                    outX[index] = dpx * scale + fx;
                    outY[index] = dpy * scale + fy;
                    outZ[index] = dpz * scale + fz;
                    outStabX[index] = stx * scale;
                    outStabY[index] = sty * scale;
                    outStabZ[index] = stz * scale;
                }
                else
                {
                    outX[index] = 0;
                    outY[index] = 0;
                    outZ[index] = 0;
                    outStabX[index] = 0;
                    outStabY[index] = 0;
                    outStabZ[index] = 0;
                }
            }
        }

        #endregion

        #region Bounds, Velocities, Sleeping

        /// <summary>
        /// Clamps particles to the bounds, derives velocities from the distance
        /// travelled, updates lifetimes and puts resting particles to sleep.
        /// </summary>
        private void Finish(float dt)
        {
            float minX = (float)Bounds.Min.X, minY = (float)Bounds.Min.Y, minZ = (float)Bounds.Min.Z;
            float maxX = (float)Bounds.Max.X, maxY = (float)Bounds.Max.Y, maxZ = (float)Bounds.Max.Z;
            float restitution = Restitution;
            float friction = Math.Max(Friction, 0);
            float invDt = 1f / dt;
            bool allowSleeping = AllowSleeping;
            float sleepSpeedSq = SleepSpeed * SleepSpeed;
            int sleepSubSteps = Math.Clamp(SleepSubSteps, 1, 255);

            // A bounce needs a real impact; resting contact must not jitter.
            const float bounceThreshold = 1.0f;

            bool solved = EnableCollisions && _maxRadius > 0 && _awakeColliders > 0;
            int active = 0;
            object gate = new();

            ForRanges(_count, (start, end) =>
            {
                int localActive = 0;

                for (int i = start; i < end; i++)
                {
                    byte flags = _flags[i];
                    if ((flags & FlagAlive) == 0)
                        continue;

                    if (_lifetime[i] != float.MaxValue)
                    {
                        _lifetime[i] -= dt;
                        if (_lifetime[i] <= 0)
                        {
                            _flags[i] = (byte)(flags & ~FlagAlive);
                            continue;
                        }
                    }

                    if ((flags & (FlagFrozen | FlagSleeping)) != 0 || _invMass[i] <= 0)
                        continue;

                    localActive++;

                    float r = _radius[i];
                    float x = _posX[i], y = _posY[i], z = _posZ[i];
                    float moveX = x - _prevX[i], moveY = y - _prevY[i], moveZ = z - _prevZ[i];
                    byte hit = 0;

                    // Bounds: push back inside and let friction eat part of the slide.
                    float depth;
                    if ((depth = minX + r - x) > 0) { x = minX + r; hit |= 1; Slide(ref y, ref z, moveY, moveZ, depth, friction); }
                    else if ((depth = x - (maxX - r)) > 0) { x = maxX - r; hit |= 1; Slide(ref y, ref z, moveY, moveZ, depth, friction); }

                    if ((depth = minY + r - y) > 0) { y = minY + r; hit |= 2; Slide(ref x, ref z, moveX, moveZ, depth, friction); }
                    else if ((depth = y - (maxY - r)) > 0) { y = maxY - r; hit |= 2; Slide(ref x, ref z, moveX, moveZ, depth, friction); }

                    if ((depth = minZ + r - z) > 0) { z = minZ + r; hit |= 4; Slide(ref x, ref y, moveX, moveY, depth, friction); }
                    else if ((depth = z - (maxZ - r)) > 0) { z = maxZ - r; hit |= 4; Slide(ref x, ref y, moveX, moveY, depth, friction); }

                    _posX[i] = x;
                    _posY[i] = y;
                    _posZ[i] = z;

                    // A particle that nothing touched keeps the velocity it was
                    // integrated with. Re-deriving it from the positions would only
                    // add rounding error, which is noticeable far from the origin.
                    bool corrected = hit != 0 || _groupOf[i] >= 0 || _softOf[i] >= 0 ||
                        (solved && (flags & FlagCollides) != 0 &&
                         (_dpX[i] != 0 || _dpY[i] != 0 || _dpZ[i] != 0));
                    if (!corrected)
                    {
                        if (allowSleeping) _rest[i] = 0;
                        continue;
                    }

                    // Velocity from the distance actually travelled.
                    float vx = (x - _prevX[i]) * invDt;
                    float vy = (y - _prevY[i]) * invDt;
                    float vz = (z - _prevZ[i]) * invDt;

                    if (solved && (flags & FlagCollides) != 0)
                    {
                        // Part of the travel was the removal of overlap that existed
                        // before this sub-step. Counting it as velocity makes packs
                        // that start out intersecting blow apart; ignoring it leaves
                        // grains in a pile with a downward velocity that is undone by
                        // the same correction every step, so the pile never rests.
                        // Of the two estimates the slower one is kept: it is right in
                        // both cases and can never add energy.
                        float ax = vx - _stabX[i] * invDt;
                        float ay = vy - _stabY[i] * invDt;
                        float az = vz - _stabZ[i] * invDt;
                        if (ax * ax + ay * ay + az * az < vx * vx + vy * vy + vz * vz)
                        {
                            vx = ax; vy = ay; vz = az;
                        }
                    }

                    if (hit != 0 && restitution > 0)
                    {
                        // _vel still holds the velocity the particle arrived with.
                        if ((hit & 1) != 0 && MathF.Abs(_velX[i]) > bounceThreshold) vx = -_velX[i] * restitution;
                        if ((hit & 2) != 0 && MathF.Abs(_velY[i]) > bounceThreshold) vy = -_velY[i] * restitution;
                        if ((hit & 4) != 0 && MathF.Abs(_velZ[i]) > bounceThreshold) vz = -_velZ[i] * restitution;
                    }

                    _velX[i] = vx;
                    _velY[i] = vy;
                    _velZ[i] = vz;

                    if (!allowSleeping || _groupOf[i] >= 0 || _softOf[i] >= 0)
                        continue;

                    if (vx * vx + vy * vy + vz * vz < sleepSpeedSq)
                    {
                        int rest = _rest[i] + 1;
                        if (rest >= sleepSubSteps)
                        {
                            _flags[i] = (byte)(flags | FlagSleeping);
                            _velX[i] = _velY[i] = _velZ[i] = 0;
                            rest = 0;
                        }
                        _rest[i] = (byte)rest;
                    }
                    else
                    {
                        _rest[i] = 0;
                    }
                }

                lock (gate) active += localActive;
            });

            _activeCount = active;
        }

        /// <summary>
        /// Coulomb friction against a wall: removes up to <c>mu * depth</c> of the
        /// tangential movement of this sub-step.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void Slide(ref float a, ref float b, float moveA, float moveB, float depth, float friction)
        {
            if (friction <= 0)
                return;

            float slideSq = moveA * moveA + moveB * moveB;
            if (slideSq < 1e-16f)
                return;

            float slide = MathF.Sqrt(slideSq);
            float limit = friction * depth;
            float scale = slide <= limit ? 1f : limit / slide;
            a -= moveA * scale;
            b -= moveB * scale;
        }

        #endregion
    }
}
