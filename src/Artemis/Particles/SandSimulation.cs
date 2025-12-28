using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using System.Numerics;
using Artemis.Core;
using Artemis.Compute;

namespace Artemis.Particles
{
    /// <summary>
    /// High-performance particle simulation for sand and granular materials.
    /// Uses multi-threading, SIMD, and optional GPU acceleration.
    /// </summary>
    public class SandSimulation
    {
        #region Fields

        // Structure of Arrays for SIMD-friendly memory layout
        private float[] _posX;
        private float[] _posY;
        private float[] _posZ;
        private float[] _velX;
        private float[] _velY;
        private float[] _velZ;
        private float[] _radii;
        private uint[] _colors;
        private byte[] _flags; // bit 0 = active, bit 1 = settled
        private int[] _groupIds;

        private int _capacity;
        private int _count;
        private int _activeCount;

        // Spatial partitioning
        private readonly Dictionary<long, List<int>> _grid;
        private readonly float _cellSize;
        private readonly float _invCellSize;

        // Thread-local storage for collision pairs
        private readonly ThreadLocal<List<(int, int)>> _threadLocalPairs;

        // GPU compute (optional)
        private GpuCompute? _gpuCompute;
        private bool _useGpu;

        // Thread synchronization
        private readonly object _countLock = new object();
        private SpinLock[] _particleLocks;

        #endregion

        #region Properties

        /// <summary>
        /// Gets all particles as a readonly view.
        /// </summary>
        public IReadOnlyList<SandParticle> Particles => GetParticleList();

        /// <summary>
        /// Gets the number of particles.
        /// </summary>
        public int ParticleCount => _count;

        /// <summary>
        /// Gets the number of active particles.
        /// </summary>
        public int ActiveParticleCount => _activeCount;

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>
        /// Gets or sets the simulation bounds.
        /// </summary>
        public AABB Bounds { get; set; }

        /// <summary>
        /// Gets or sets the friction coefficient.
        /// </summary>
        public double Friction { get; set; } = 0.3;

        /// <summary>
        /// Gets or sets the restitution (bounciness).
        /// </summary>
        public double Restitution { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets the default particle radius.
        /// </summary>
        public double ParticleRadius { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets whether to use GPU acceleration.
        /// </summary>
        public bool UseGPU
        {
            get => _useGpu;
            set
            {
                if (value && _gpuCompute == null)
                {
                    InitializeGpu();
                }
                _useGpu = value && _gpuCompute != null;
            }
        }

        /// <summary>
        /// Gets or sets the number of threads (0 = auto).
        /// </summary>
        public int ThreadCount { get; set; } = 0;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new sand simulation.
        /// </summary>
        public SandSimulation(AABB bounds, double cellSize = 0.5)
        {
            Bounds = bounds;
            _cellSize = (float)cellSize;
            _invCellSize = 1f / _cellSize;

            _capacity = 1024;
            _count = 0;
            _activeCount = 0;

            // Allocate SoA arrays
            AllocateArrays(_capacity);

            _grid = new Dictionary<long, List<int>>();
            _threadLocalPairs = new ThreadLocal<List<(int, int)>>(() => new List<(int, int)>(256), true);
            _particleLocks = new SpinLock[_capacity];
            for (int i = 0; i < _capacity; i++)
                _particleLocks[i] = new SpinLock(false);
        }

        private void AllocateArrays(int capacity)
        {
            _posX = new float[capacity];
            _posY = new float[capacity];
            _posZ = new float[capacity];
            _velX = new float[capacity];
            _velY = new float[capacity];
            _velZ = new float[capacity];
            _radii = new float[capacity];
            _colors = new uint[capacity];
            _flags = new byte[capacity];
            _groupIds = new int[capacity];
        }

        private void EnsureCapacity(int required)
        {
            if (required <= _capacity) return;

            int newCapacity = Math.Max(required, _capacity * 2);

            var newPosX = new float[newCapacity];
            var newPosY = new float[newCapacity];
            var newPosZ = new float[newCapacity];
            var newVelX = new float[newCapacity];
            var newVelY = new float[newCapacity];
            var newVelZ = new float[newCapacity];
            var newRadii = new float[newCapacity];
            var newColors = new uint[newCapacity];
            var newFlags = new byte[newCapacity];
            var newGroupIds = new int[newCapacity];

            Array.Copy(_posX, newPosX, _count);
            Array.Copy(_posY, newPosY, _count);
            Array.Copy(_posZ, newPosZ, _count);
            Array.Copy(_velX, newVelX, _count);
            Array.Copy(_velY, newVelY, _count);
            Array.Copy(_velZ, newVelZ, _count);
            Array.Copy(_radii, newRadii, _count);
            Array.Copy(_colors, newColors, _count);
            Array.Copy(_flags, newFlags, _count);
            Array.Copy(_groupIds, newGroupIds, _count);

            _posX = newPosX;
            _posY = newPosY;
            _posZ = newPosZ;
            _velX = newVelX;
            _velY = newVelY;
            _velZ = newVelZ;
            _radii = newRadii;
            _colors = newColors;
            _flags = newFlags;
            _groupIds = newGroupIds;

            var newLocks = new SpinLock[newCapacity];
            for (int i = 0; i < newCapacity; i++)
                newLocks[i] = new SpinLock(false);
            _particleLocks = newLocks;

            _capacity = newCapacity;
        }

        #endregion

        #region GPU Initialization

        private void InitializeGpu()
        {
            try
            {
                _gpuCompute = new GpuCompute();
                _gpuCompute.Initialize();
            }
            catch
            {
                _gpuCompute = null;
            }
        }

        #endregion

        #region Particle Management

        /// <summary>
        /// Adds a sand particle.
        /// </summary>
        public int AddParticle(Vector3D position, uint color, double? radius = null)
        {
            lock (_countLock)
            {
                EnsureCapacity(_count + 1);

                int index = _count;
                _posX[index] = (float)position.X;
                _posY[index] = (float)position.Y;
                _posZ[index] = (float)position.Z;
                _velX[index] = 0;
                _velY[index] = 0;
                _velZ[index] = 0;
                _radii[index] = (float)(radius ?? ParticleRadius);
                _colors[index] = color;
                _flags[index] = 1; // Active, not settled
                _groupIds[index] = -1;

                _count++;
                _activeCount++;
                return index;
            }
        }

        /// <summary>
        /// Adds a spherical cluster of sand particles.
        /// </summary>
        public (int startIndex, int count) AddSandBall(
            Vector3D center,
            double ballRadius,
            uint color,
            double? particleRadius = null)
        {
            float pRadius = (float)(particleRadius ?? ParticleRadius);
            float spacing = pRadius * 2.1f;

            int startIndex = _count;
            int steps = (int)Math.Ceiling(ballRadius / spacing);
            var positions = new List<(float x, float y, float z)>();

            // Pre-calculate positions
            for (int x = -steps; x <= steps; x++)
            {
                for (int y = -steps; y <= steps; y++)
                {
                    for (int z = -steps; z <= steps; z++)
                    {
                        float ox = x * spacing;
                        float oy = y * spacing;
                        float oz = z * spacing;
                        float dist = MathF.Sqrt(ox * ox + oy * oy + oz * oz);
                        if (dist <= (float)ballRadius - pRadius)
                        {
                            positions.Add(((float)center.X + ox, (float)center.Y + oy, (float)center.Z + oz));
                        }
                    }
                }
            }

            lock (_countLock)
            {
                EnsureCapacity(_count + positions.Count);

                foreach (var (px, py, pz) in positions)
                {
                    int index = _count;
                    _posX[index] = px;
                    _posY[index] = py;
                    _posZ[index] = pz;
                    _velX[index] = 0;
                    _velY[index] = 0;
                    _velZ[index] = 0;
                    _radii[index] = pRadius;
                    _colors[index] = color;
                    _flags[index] = 1;
                    _groupIds[index] = -1;
                    _count++;
                    _activeCount++;
                }
            }

            return (startIndex, positions.Count);
        }

        /// <summary>
        /// Removes a particle by index.
        /// </summary>
        public void RemoveParticle(int index)
        {
            if (index >= 0 && index < _count && (_flags[index] & 1) != 0)
            {
                _flags[index] = 0;
                Interlocked.Decrement(ref _activeCount);
            }
        }

        /// <summary>
        /// Removes all particles with the specified group ID.
        /// </summary>
        public int RemoveGroup(int groupId)
        {
            int removed = 0;
            for (int i = 0; i < _count; i++)
            {
                if (_groupIds[i] == groupId && (_flags[i] & 1) != 0)
                {
                    _flags[i] = 0;
                    removed++;
                }
            }
            Interlocked.Add(ref _activeCount, -removed);
            return removed;
        }

        /// <summary>
        /// Clears all particles.
        /// </summary>
        public void Clear()
        {
            lock (_countLock)
            {
                _count = 0;
                _activeCount = 0;
                _grid.Clear();
            }
        }

        private List<SandParticle> GetParticleList()
        {
            var result = new List<SandParticle>(_count);
            for (int i = 0; i < _count; i++)
            {
                result.Add(new SandParticle
                {
                    Position = new Vector3D(_posX[i], _posY[i], _posZ[i]),
                    PreviousPosition = new Vector3D(_posX[i], _posY[i], _posZ[i]),
                    Velocity = new Vector3D(_velX[i], _velY[i], _velZ[i]),
                    Radius = _radii[i],
                    Color = _colors[i],
                    IsActive = (_flags[i] & 1) != 0,
                    IsSettled = (_flags[i] & 2) != 0,
                    GroupId = _groupIds[i]
                });
            }
            return result;
        }

        #endregion

        #region Simulation

        /// <summary>
        /// Updates the simulation with multi-threading and SIMD.
        /// </summary>
        public void Update(double deltaTime, int subSteps = 4)
        {
            float dt = (float)deltaTime / subSteps;
            float gravX = (float)Gravity.X;
            float gravY = (float)Gravity.Y;
            float gravZ = (float)Gravity.Z;

            for (int step = 0; step < subSteps; step++)
            {
                // 1. Rebuild spatial grid (parallel)
                RebuildGridParallel();

                // 2. Apply gravity and integrate (SIMD + parallel)
                IntegrateParallel(dt, gravX, gravY, gravZ);

                // 3. Handle collisions (parallel)
                HandleCollisionsParallel();

                // 4. Handle bounds (parallel + SIMD)
                HandleBoundsParallel();

                // 5. Check settling (parallel)
                CheckSettlingParallel();
            }
        }

        #endregion

        #region Parallel Grid Rebuild

        private void RebuildGridParallel()
        {
            _grid.Clear();

            // Single-threaded grid rebuild (dictionary not thread-safe for writes)
            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & 1) == 0) continue;

                long cellKey = GetCellKey(_posX[i], _posY[i], _posZ[i]);
                if (!_grid.TryGetValue(cellKey, out var list))
                {
                    list = new List<int>(16);
                    _grid[cellKey] = list;
                }
                list.Add(i);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long GetCellKey(float x, float y, float z)
        {
            int cx = (int)MathF.Floor(x * _invCellSize);
            int cy = (int)MathF.Floor(y * _invCellSize);
            int cz = (int)MathF.Floor(z * _invCellSize);
            // Pack into 64 bits: 21 bits each for x, y, z
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (cz & 0x1FFFFF);
        }

        #endregion

        #region SIMD Integration

        private void IntegrateParallel(float dt, float gravX, float gravY, float gravZ)
        {
            int simdWidth = Vector<float>.Count;
            var dtVec = new Vector<float>(dt);
            var gravXVec = new Vector<float>(gravX * dt);
            var gravYVec = new Vector<float>(gravY * dt);
            var gravZVec = new Vector<float>(gravZ * dt);

            // Process in SIMD-width chunks
            int chunks = (_count + simdWidth - 1) / simdWidth;

            Parallel.For(0, chunks, chunk =>
            {
                int start = chunk * simdWidth;
                int end = Math.Min(start + simdWidth, _count);

                if (end - start == simdWidth && start + simdWidth <= _count)
                {
                    // Full SIMD path
                    var velX = new Vector<float>(_velX, start);
                    var velY = new Vector<float>(_velY, start);
                    var velZ = new Vector<float>(_velZ, start);
                    var posX = new Vector<float>(_posX, start);
                    var posY = new Vector<float>(_posY, start);
                    var posZ = new Vector<float>(_posZ, start);
                    var flags = new Vector<byte>(_flags, start);

                    // Apply gravity to velocity
                    velX += gravXVec;
                    velY += gravYVec;
                    velZ += gravZVec;

                    // Integrate position
                    posX += velX * dtVec;
                    posY += velY * dtVec;
                    posZ += velZ * dtVec;

                    // Store back
                    velX.CopyTo(_velX, start);
                    velY.CopyTo(_velY, start);
                    velZ.CopyTo(_velZ, start);
                    posX.CopyTo(_posX, start);
                    posY.CopyTo(_posY, start);
                    posZ.CopyTo(_posZ, start);
                }
                else
                {
                    // Scalar fallback for remainder
                    for (int i = start; i < end; i++)
                    {
                        if ((_flags[i] & 1) == 0 || (_flags[i] & 2) != 0) continue;

                        _velX[i] += gravX * dt;
                        _velY[i] += gravY * dt;
                        _velZ[i] += gravZ * dt;

                        _posX[i] += _velX[i] * dt;
                        _posY[i] += _velY[i] * dt;
                        _posZ[i] += _velZ[i] * dt;
                    }
                }
            });
        }

        #endregion

        #region Parallel Collision Detection

        private void HandleCollisionsParallel()
        {
            float restitution = (float)Restitution;
            float friction = (float)Friction;

            // Collect collision pairs from all cells in parallel
            var cells = _grid.ToArray();

            Parallel.ForEach(cells, cell =>
            {
                var indices = cell.Value;
                var cellKey = cell.Key;

                // Unpack cell coordinates
                int cx = (int)((cellKey >> 42) & 0x1FFFFF);
                int cy = (int)((cellKey >> 21) & 0x1FFFFF);
                int cz = (int)(cellKey & 0x1FFFFF);
                if (cx >= 0x100000) cx -= 0x200000; // Sign extend
                if (cy >= 0x100000) cy -= 0x200000;
                if (cz >= 0x100000) cz -= 0x200000;

                // Check within cell
                for (int i = 0; i < indices.Count; i++)
                {
                    for (int j = i + 1; j < indices.Count; j++)
                    {
                        ResolveCollisionSafe(indices[i], indices[j], restitution, friction);
                    }
                }

                // Check neighboring cells (only positive direction to avoid duplicates)
                for (int dx = 0; dx <= 1; dx++)
                {
                    for (int dy = 0; dy <= 1; dy++)
                    {
                        for (int dz = 0; dz <= 1; dz++)
                        {
                            if (dx == 0 && dy == 0 && dz == 0) continue;

                            long neighborKey = GetCellKeyFromCoords(cx + dx, cy + dy, cz + dz);
                            if (!_grid.TryGetValue(neighborKey, out var neighborIndices)) continue;

                            foreach (var i in indices)
                            {
                                foreach (var j in neighborIndices)
                                {
                                    ResolveCollisionSafe(i, j, restitution, friction);
                                }
                            }
                        }
                    }
                }
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long GetCellKeyFromCoords(int cx, int cy, int cz)
        {
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (cz & 0x1FFFFF);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void ResolveCollisionSafe(int i, int j, float restitution, float friction)
        {
            // Ensure consistent lock ordering to prevent deadlocks
            int first = Math.Min(i, j);
            int second = Math.Max(i, j);

            bool lockTakenFirst = false;
            bool lockTakenSecond = false;

            try
            {
                _particleLocks[first].Enter(ref lockTakenFirst);
                _particleLocks[second].Enter(ref lockTakenSecond);

                if ((_flags[i] & 1) == 0 || (_flags[j] & 1) == 0) return;

                float dx = _posX[j] - _posX[i];
                float dy = _posY[j] - _posY[i];
                float dz = _posZ[j] - _posZ[i];

                float distSq = dx * dx + dy * dy + dz * dz;
                float minDist = _radii[i] + _radii[j];

                if (distSq >= minDist * minDist || distSq < 1e-8f) return;

                float dist = MathF.Sqrt(distSq);
                float invDist = 1f / dist;
                float nx = dx * invDist;
                float ny = dy * invDist;
                float nz = dz * invDist;
                float penetration = minDist - dist;

                // Position correction
                float correction = penetration * 0.5f;
                _posX[i] -= nx * correction;
                _posY[i] -= ny * correction;
                _posZ[i] -= nz * correction;
                _posX[j] += nx * correction;
                _posY[j] += ny * correction;
                _posZ[j] += nz * correction;

                // Velocity correction
                float dvx = _velX[j] - _velX[i];
                float dvy = _velY[j] - _velY[i];
                float dvz = _velZ[j] - _velZ[i];

                float velAlongNormal = dvx * nx + dvy * ny + dvz * nz;
                if (velAlongNormal > 0) return;

                float impulseMag = -(1 + restitution) * velAlongNormal * 0.5f;
                float impulseX = nx * impulseMag;
                float impulseY = ny * impulseMag;
                float impulseZ = nz * impulseMag;

                _velX[i] -= impulseX;
                _velY[i] -= impulseY;
                _velZ[i] -= impulseZ;
                _velX[j] += impulseX;
                _velY[j] += impulseY;
                _velZ[j] += impulseZ;

                // Tangential friction
                float tanX = dvx - nx * velAlongNormal;
                float tanY = dvy - ny * velAlongNormal;
                float tanZ = dvz - nz * velAlongNormal;
                float tanLenSq = tanX * tanX + tanY * tanY + tanZ * tanZ;

                if (tanLenSq > 1e-8f)
                {
                    float tanLen = MathF.Sqrt(tanLenSq);
                    float invTanLen = 1f / tanLen;
                    tanX *= invTanLen;
                    tanY *= invTanLen;
                    tanZ *= invTanLen;

                    float frictionImpulse = friction * impulseMag;
                    _velX[i] += tanX * frictionImpulse;
                    _velY[i] += tanY * frictionImpulse;
                    _velZ[i] += tanZ * frictionImpulse;
                    _velX[j] -= tanX * frictionImpulse;
                    _velY[j] -= tanY * frictionImpulse;
                    _velZ[j] -= tanZ * frictionImpulse;
                }

                // Wake up settled particles
                _flags[i] &= unchecked((byte)~2);
                _flags[j] &= unchecked((byte)~2);
            }
            finally
            {
                if (lockTakenSecond) _particleLocks[second].Exit(false);
                if (lockTakenFirst) _particleLocks[first].Exit(false);
            }
        }

        #endregion

        #region Parallel Bounds Handling

        private void HandleBoundsParallel()
        {
            float minX = (float)Bounds.Min.X;
            float minY = (float)Bounds.Min.Y;
            float minZ = (float)Bounds.Min.Z;
            float maxX = (float)Bounds.Max.X;
            float maxY = (float)Bounds.Max.Y;
            float maxZ = (float)Bounds.Max.Z;
            float restitution = (float)Restitution;
            float friction = 1f - (float)Friction;

            Parallel.For(0, _count, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                float r = _radii[i];

                // Floor
                if (_posY[i] - r < minY)
                {
                    _posY[i] = minY + r;
                    _velY[i] = -_velY[i] * restitution;
                    _velX[i] *= friction;
                    _velZ[i] *= friction;
                }
                // Ceiling
                if (_posY[i] + r > maxY)
                {
                    _posY[i] = maxY - r;
                    _velY[i] = -_velY[i] * restitution;
                }
                // Walls X
                if (_posX[i] - r < minX)
                {
                    _posX[i] = minX + r;
                    _velX[i] = -_velX[i] * restitution;
                }
                if (_posX[i] + r > maxX)
                {
                    _posX[i] = maxX - r;
                    _velX[i] = -_velX[i] * restitution;
                }
                // Walls Z
                if (_posZ[i] - r < minZ)
                {
                    _posZ[i] = minZ + r;
                    _velZ[i] = -_velZ[i] * restitution;
                }
                if (_posZ[i] + r > maxZ)
                {
                    _posZ[i] = maxZ - r;
                    _velZ[i] = -_velZ[i] * restitution;
                }
            });
        }

        #endregion

        #region Parallel Settling

        private void CheckSettlingParallel()
        {
            const float settleThresholdSq = 0.0001f; // 0.01^2

            Parallel.For(0, _count, i =>
            {
                if ((_flags[i] & 1) == 0 || (_flags[i] & 2) != 0) return;

                float velSq = _velX[i] * _velX[i] + _velY[i] * _velY[i] + _velZ[i] * _velZ[i];
                if (velSq < settleThresholdSq)
                {
                    _velX[i] = 0;
                    _velY[i] = 0;
                    _velZ[i] = 0;
                    _flags[i] |= 2; // Mark as settled
                }
            });
        }

        #endregion

        #region Queries

        /// <summary>
        /// Finds connected particles of the same color using flood fill.
        /// </summary>
        public List<int> FindConnectedParticles(int startIndex)
        {
            if (startIndex < 0 || startIndex >= _count || (_flags[startIndex] & 1) == 0)
                return new List<int>();

            var result = new List<int>();
            var visited = new HashSet<int>();
            var queue = new Queue<int>();
            uint targetColor = _colors[startIndex];

            queue.Enqueue(startIndex);
            visited.Add(startIndex);

            float connectionDistance = (float)ParticleRadius * 2.5f;
            float connectionDistSq = connectionDistance * connectionDistance;

            while (queue.Count > 0)
            {
                int current = queue.Dequeue();
                result.Add(current);

                float cx = _posX[current];
                float cy = _posY[current];
                float cz = _posZ[current];

                // Use spatial hash for faster neighbor lookup
                int cellX = (int)MathF.Floor(cx * _invCellSize);
                int cellY = (int)MathF.Floor(cy * _invCellSize);
                int cellZ = (int)MathF.Floor(cz * _invCellSize);

                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            long cellKey = GetCellKeyFromCoords(cellX + dx, cellY + dy, cellZ + dz);
                            if (!_grid.TryGetValue(cellKey, out var indices)) continue;

                            foreach (int i in indices)
                            {
                                if (visited.Contains(i) || (_flags[i] & 1) == 0) continue;
                                if (_colors[i] != targetColor) continue;

                                float ddx = _posX[i] - cx;
                                float ddy = _posY[i] - cy;
                                float ddz = _posZ[i] - cz;
                                float distSq = ddx * ddx + ddy * ddy + ddz * ddz;

                                if (distSq <= connectionDistSq)
                                {
                                    visited.Add(i);
                                    queue.Enqueue(i);
                                }
                            }
                        }
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Checks if particles span from one side to another.
        /// </summary>
        public bool CheckSpansAxis(uint color, int axis)
        {
            float minThreshold = GetAxisMin(axis) + (float)ParticleRadius * 3;
            float maxThreshold = GetAxisMax(axis) - (float)ParticleRadius * 3;

            var startParticles = new List<int>();
            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & 1) == 0 || _colors[i] != color) continue;

                float pos = GetAxisValue(i, axis);
                if (pos <= minThreshold)
                    startParticles.Add(i);
            }

            foreach (int start in startParticles)
            {
                var connected = FindConnectedParticles(start);
                foreach (int idx in connected)
                {
                    float pos = GetAxisValue(idx, axis);
                    if (pos >= maxThreshold)
                        return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Checks if there's a hole at the bottom.
        /// </summary>
        public bool HasHoleAtBottom()
        {
            float bottom = (float)Bounds.Min.Y + (float)ParticleRadius * 2;
            float checkRadius = (float)ParticleRadius * 2;
            float checkRadiusSq = checkRadius * checkRadius;
            float stepX = (float)ParticleRadius * 2;
            float stepZ = (float)ParticleRadius * 2;

            for (float x = (float)Bounds.Min.X + checkRadius; x < (float)Bounds.Max.X - checkRadius; x += stepX)
            {
                for (float z = (float)Bounds.Min.Z + checkRadius; z < (float)Bounds.Max.Z - checkRadius; z += stepZ)
                {
                    bool hasParticle = false;
                    for (int i = 0; i < _count && !hasParticle; i++)
                    {
                        if ((_flags[i] & 1) == 0) continue;

                        float dx = _posX[i] - x;
                        float dy = _posY[i] - bottom;
                        float dz = _posZ[i] - z;
                        if (dx * dx + dy * dy + dz * dz < checkRadiusSq)
                            hasParticle = true;
                    }

                    if (!hasParticle)
                        return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Checks if the terrarium is full.
        /// </summary>
        public bool IsFull()
        {
            float topThreshold = (float)Bounds.Max.Y - (float)ParticleRadius * 4;

            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & 1) != 0 && _posY[i] >= topThreshold)
                    return true;
            }

            return false;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private float GetAxisValue(int index, int axis) => axis switch
        {
            0 => _posX[index],
            1 => _posY[index],
            2 => _posZ[index],
            _ => 0
        };

        private float GetAxisMin(int axis) => axis switch
        {
            0 => (float)Bounds.Min.X,
            1 => (float)Bounds.Min.Y,
            2 => (float)Bounds.Min.Z,
            _ => 0
        };

        private float GetAxisMax(int axis) => axis switch
        {
            0 => (float)Bounds.Max.X,
            1 => (float)Bounds.Max.Y,
            2 => (float)Bounds.Max.Z,
            _ => 0
        };

        #endregion
    }

    /// <summary>
    /// Represents a sand particle.
    /// </summary>
    public struct SandParticle
    {
        public Vector3D Position;
        public Vector3D PreviousPosition;
        public Vector3D Velocity;
        public double Radius;
        public uint Color;
        public bool IsActive;
        public bool IsSettled;
        public int GroupId;
    }
}
