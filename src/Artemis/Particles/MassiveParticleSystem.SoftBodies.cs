using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;
using Artemis.Compute;
using Quaternion = System.Numerics.Quaternion;

namespace Artemis.Particles
{
    // Soft bodies: jelly, rubber, flesh.
    //
    // A soft body is a set of particles covered by overlapping clusters (clustered
    // shape matching, Mueller et al. 2005, as used by NVIDIA FleX). Every sub-step,
    // each cluster finds the rigid placement of its rest shape that best fits where
    // its particles are, and each particle is pulled part of the way towards the
    // mean of its clusters' goals. The fraction it is pulled is the stiffness. Small
    // clusters let the body bend and wobble locally; because every goal comes from a
    // rotated rest shape, the body keeps its volume and always springs back.
    public partial class MassiveParticleSystem
    {
        #region Types and Fields

        private sealed class SoftBody
        {
            public int[] Members = Array.Empty<int>();
            public float Stiffness;
            public float Damping;
            public bool Alive = true;
            public List<Cluster> Clusters = new();
        }

        private sealed class Cluster
        {
            public int[] Members = Array.Empty<int>();
            public float[] RestX = Array.Empty<float>();
            public float[] RestY = Array.Empty<float>();
            public float[] RestZ = Array.Empty<float>();
            public Quaternion Rotation = Quaternion.Identity;
        }

        private readonly List<SoftBody> _softBodies = new();
        private int[] _softOf = Array.Empty<int>();
        private int _aliveSoftBodies;

        // Clusters of all soft bodies, flattened (rebuilt when the structure changes).
        private bool _softStructureDirty = true;
        private int _softBuiltForCount = -1;
        private int _clusterCount;
        private int _slotCount;
        private int[] _clusterStart = Array.Empty<int>();
        private int[] _clusterMembers = Array.Empty<int>();
        private float[] _clusterRestX = Array.Empty<float>();
        private float[] _clusterRestY = Array.Empty<float>();
        private float[] _clusterRestZ = Array.Empty<float>();
        private float[] _clusterRotation = Array.Empty<float>();  // x, y, z, w per cluster
        private int[] _slotStart = Array.Empty<int>();             // per particle, into _slots
        private int[] _slots = Array.Empty<int>();
        private float[] _goalX = Array.Empty<float>();
        private float[] _goalY = Array.Empty<float>();
        private float[] _goalZ = Array.Empty<float>();
        private float[] _softAlpha = Array.Empty<float>();
        private float[] _softBeta = Array.Empty<float>();

        // Body-level damping: member lists in chunks, and the rigid motion of each body.
        private const int DampChunk = 256;
        private int[] _bodyStart = Array.Empty<int>();      // per body, into _bodyMembers
        private int[] _bodyMembers = Array.Empty<int>();
        private int[] _bodyChunkStart = Array.Empty<int>(); // per body, into the chunk arrays
        private int[] _chunkBody = Array.Empty<int>();
        private int[] _chunkFrom = Array.Empty<int>();
        private int[] _chunkTo = Array.Empty<int>();
        private int _chunkCount;
        private float[] _chunkPartial = Array.Empty<float>(); // 16 sums per chunk
        private float[] _bodyRigid = Array.Empty<float>();    // vcm, centre, omega per body
        private bool _softUploaded;

        #endregion

        #region Properties

        /// <summary>
        /// Number of soft bodies.
        /// </summary>
        public int SoftBodyCount => _aliveSoftBodies;

        #endregion

        #region Creation

        private void AllocateSoftArrays(int capacity)
        {
            _softOf = new int[capacity];
            for (int i = 0; i < capacity; i++)
                _softOf[i] = -1;
        }

        /// <summary>
        /// Rebuilds a mesh out of particles that form one soft body.
        /// </summary>
        /// <param name="mesh">The mesh, already placed where it should appear. It should be closed.</param>
        /// <param name="particleRadius">Radius of the particles.</param>
        /// <param name="softBodyId">Id of the soft body, or -1 if nothing was created.</param>
        /// <param name="stiffness">
        /// 0 to 1: how firmly the body returns to its shape. About 0.05 is jelly,
        /// 0.3 rubber, 0.9 and above nearly rigid.
        /// </param>
        /// <param name="clusterRadius">
        /// Size of the clusters in particle diameters. Smaller clusters bend more
        /// locally; larger ones keep the overall shape better. 1.5 to 3 is the useful range.
        /// </param>
        /// <param name="mass">Mass of one particle.</param>
        /// <param name="color">Colour of the particles.</param>
        /// <param name="damping">
        /// 0 to 1: internal friction. It damps deformation - squashing, stretching,
        /// wobbling and the rebound they cause - but not the motion of the body as a
        /// whole (falling, sliding, spinning). 0 bounces like a rubber ball; around
        /// 0.5 lands with hardly any bounce.
        /// </param>
        /// <returns>The number of particles created.</returns>
        public int SpawnSoftBodyFromMesh(
            TriangleMesh mesh,
            float particleRadius,
            out int softBodyId,
            float stiffness = 0.3f,
            float clusterRadius = 2f,
            float mass = 1.0f,
            uint color = 0xFFFFFFFF,
            float damping = 0.1f)
        {
            int created = SpawnFromMesh(mesh, particleRadius, out int first, solid: true, frozen: false, mass, color, jitter: 0);
            softBodyId = created > 0 ? CreateSoftBody(first, created, stiffness, clusterRadius, damping) : -1;
            return created;
        }

        /// <summary>
        /// Turns a range of particles into one soft body. The current positions are
        /// the rest shape.
        /// </summary>
        /// <param name="firstIndex">First particle of the range.</param>
        /// <param name="count">Number of particles.</param>
        /// <param name="stiffness">See <see cref="SpawnSoftBodyFromMesh"/>.</param>
        /// <param name="clusterRadius">See <see cref="SpawnSoftBodyFromMesh"/>.</param>
        /// <param name="damping">See <see cref="SpawnSoftBodyFromMesh"/>.</param>
        /// <returns>The id of the soft body.</returns>
        public int CreateSoftBody(int firstIndex, int count, float stiffness = 0.3f, float clusterRadius = 2f, float damping = 0.1f)
        {
            TouchHost();
            if (firstIndex < 0 || count <= 0 || firstIndex + count > _count)
                throw new ArgumentOutOfRangeException(nameof(count));
            if (!(clusterRadius > 0))
                throw new ArgumentOutOfRangeException(nameof(clusterRadius));

            var members = new int[count];
            float radiusSum = 0;
            for (int k = 0; k < count; k++)
            {
                int i = firstIndex + k;
                if (_groupOf[i] >= 0 || _softOf[i] >= 0)
                    throw new InvalidOperationException($"Particle {i} already belongs to another body.");
                members[k] = i;
                radiusSum += _radius[i];
            }

            var body = new SoftBody { Members = members, Stiffness = Math.Clamp(stiffness, 0f, 1f), Damping = Math.Clamp(damping, 0f, 1f) };
            int id = _softBodies.Count;
            _softBodies.Add(body);
            _aliveSoftBodies++;

            foreach (int i in members)
            {
                _softOf[i] = id;
                _flags[i] &= unchecked((byte)~(FlagFrozen | FlagSleeping));
                _rest[i] = 0;
            }

            float diameter = 2f * radiusSum / count;
            BuildClusters(body, clusterRadius * diameter);
            _softStructureDirty = true;
            return id;
        }

        /// <summary>
        /// Covers the body with overlapping spherical clusters centred on a grid whose
        /// spacing equals the cluster radius, so every particle is in at least one.
        /// </summary>
        private void BuildClusters(SoftBody body, float clusterRadius)
        {
            var members = body.Members;
            float radiusSq = clusterRadius * clusterRadius;
            float inv = 1f / clusterRadius;

            // Bucket the members on a grid of the cluster radius.
            var grid = new Dictionary<(int, int, int), List<int>>();
            var min = new Vector3(float.MaxValue);
            var max = new Vector3(float.MinValue);
            foreach (int i in members)
            {
                var p = new Vector3(_posX[i], _posY[i], _posZ[i]);
                min = Vector3.Min(min, p);
                max = Vector3.Max(max, p);
                var key = ((int)MathF.Floor(p.X * inv), (int)MathF.Floor(p.Y * inv), (int)MathF.Floor(p.Z * inv));
                if (!grid.TryGetValue(key, out var list))
                    grid[key] = list = new List<int>();
                list.Add(i);
            }

            var covered = new HashSet<int>();
            var clusterLists = new List<List<int>>();
            // Centre the grid of cluster centres on the body: a grid anchored at one
            // corner would make the body respond differently on each side, and it
            // would drift or tilt every time it is squeezed.
            int nx = (int)MathF.Ceiling((max.X - min.X) * inv) + 1;
            int ny = (int)MathF.Ceiling((max.Y - min.Y) * inv) + 1;
            int nz = (int)MathF.Ceiling((max.Z - min.Z) * inv) + 1;
            var origin = (min + max) * 0.5f - new Vector3(nx - 1, ny - 1, nz - 1) * (clusterRadius * 0.5f);

            for (int gz = 0; gz < nz; gz++)
            {
                for (int gy = 0; gy < ny; gy++)
                {
                    for (int gx = 0; gx < nx; gx++)
                    {
                        var center = origin + new Vector3(gx, gy, gz) * clusterRadius;
                        var list = new List<int>();
                        int cx = (int)MathF.Floor(center.X * inv), cy = (int)MathF.Floor(center.Y * inv), cz = (int)MathF.Floor(center.Z * inv);
                        for (int oz = -1; oz <= 1; oz++)
                            for (int oy = -1; oy <= 1; oy++)
                                for (int ox = -1; ox <= 1; ox++)
                                {
                                    if (!grid.TryGetValue((cx + ox, cy + oy, cz + oz), out var cell))
                                        continue;
                                    foreach (int i in cell)
                                    {
                                        if (Vector3.DistanceSquared(new Vector3(_posX[i], _posY[i], _posZ[i]), center) <= radiusSq)
                                            list.Add(i);
                                    }
                                }

                        // Fewer than four points cannot define a rotation reliably.
                        if (list.Count < 4)
                            continue;

                        list.Sort();
                        clusterLists.Add(list);
                        foreach (int i in list)
                            covered.Add(i);
                    }
                }
            }

            // Particles no cluster reached (thin features) join the nearest cluster.
            if (clusterLists.Count == 0)
            {
                clusterLists.Add(new List<int>(members));
                foreach (int i in members) covered.Add(i);
            }

            foreach (int i in members)
            {
                if (covered.Contains(i))
                    continue;

                var p = new Vector3(_posX[i], _posY[i], _posZ[i]);
                int best = 0;
                float bestDistance = float.MaxValue;
                for (int c = 0; c < clusterLists.Count; c++)
                {
                    var first = clusterLists[c][0];
                    float d = Vector3.DistanceSquared(p, new Vector3(_posX[first], _posY[first], _posZ[first]));
                    if (d < bestDistance) { bestDistance = d; best = c; }
                }
                clusterLists[best].Add(i);
            }

            body.Clusters.Clear();
            foreach (var list in clusterLists)
                body.Clusters.Add(MakeCluster(list.ToArray()));
        }

        /// <summary>
        /// Rest shape of a cluster: member offsets from its mass-weighted centre.
        /// </summary>
        private Cluster MakeCluster(int[] members)
        {
            double cx = 0, cy = 0, cz = 0, total = 0;
            foreach (int i in members)
            {
                double m = _mass[i];
                cx += m * _posX[i]; cy += m * _posY[i]; cz += m * _posZ[i];
                total += m;
            }
            if (total <= 0) total = 1;

            var cluster = new Cluster
            {
                Members = members,
                RestX = new float[members.Length],
                RestY = new float[members.Length],
                RestZ = new float[members.Length]
            };

            for (int k = 0; k < members.Length; k++)
            {
                int i = members[k];
                cluster.RestX[k] = _posX[i] - (float)(cx / total);
                cluster.RestY[k] = _posY[i] - (float)(cy / total);
                cluster.RestZ[k] = _posZ[i] - (float)(cz / total);
            }

            return cluster;
        }

        #endregion

        #region Queries and Control

        /// <summary>Whether a soft body exists.</summary>
        public bool IsSoftBodyAlive(int softBodyId)
            => softBodyId >= 0 && softBodyId < _softBodies.Count && _softBodies[softBodyId].Alive;

        /// <summary>Id of the soft body a particle belongs to, or -1.</summary>
        public int GetParticleSoftBody(int index) => _softOf[index];

        /// <summary>Number of particles in a soft body.</summary>
        public int GetSoftBodyParticleCount(int softBodyId) => GetSoftBody(softBodyId).Members.Length;

        /// <summary>Number of clusters a soft body is made of.</summary>
        public int GetSoftBodyClusterCount(int softBodyId) => GetSoftBody(softBodyId).Clusters.Count;

        /// <summary>Stiffness of a soft body (0 to 1).</summary>
        public float GetSoftBodyStiffness(int softBodyId) => GetSoftBody(softBodyId).Stiffness;

        /// <summary>
        /// Changes the stiffness of a soft body (0 to 1). Can be animated, for example
        /// to make a body go limp.
        /// </summary>
        public void SetSoftBodyStiffness(int softBodyId, float stiffness)
            => GetSoftBody(softBodyId).Stiffness = Math.Clamp(stiffness, 0f, 1f);

        /// <summary>Internal damping of a soft body (0 to 1).</summary>
        public float GetSoftBodyDamping(int softBodyId) => GetSoftBody(softBodyId).Damping;

        /// <summary>Changes the internal damping of a soft body (0 to 1).</summary>
        public void SetSoftBodyDamping(int softBodyId, float damping)
            => GetSoftBody(softBodyId).Damping = Math.Clamp(damping, 0f, 1f);

        /// <summary>
        /// Mass-weighted centre of a soft body.
        /// </summary>
        public Vector3 GetSoftBodyCenter(int softBodyId)
        {
            SyncHost();
            var body = GetSoftBody(softBodyId);
            double x = 0, y = 0, z = 0, total = 0;
            foreach (int i in body.Members)
            {
                double m = _mass[i];
                x += m * _posX[i]; y += m * _posY[i]; z += m * _posZ[i];
                total += m;
            }
            return total > 0 ? new Vector3((float)(x / total), (float)(y / total), (float)(z / total)) : Vector3.Zero;
        }

        /// <summary>
        /// Mass-weighted velocity of a soft body.
        /// </summary>
        public Vector3 GetSoftBodyVelocity(int softBodyId)
        {
            SyncHost();
            var body = GetSoftBody(softBodyId);
            double x = 0, y = 0, z = 0, total = 0;
            foreach (int i in body.Members)
            {
                double m = _mass[i];
                x += m * _velX[i]; y += m * _velY[i]; z += m * _velZ[i];
                total += m;
            }
            return total > 0 ? new Vector3((float)(x / total), (float)(y / total), (float)(z / total)) : Vector3.Zero;
        }

        /// <summary>
        /// Sets the velocity of a soft body: a linear velocity and, optionally, a spin
        /// (angular velocity in rad/s about its centre).
        /// </summary>
        public void SetSoftBodyVelocity(int softBodyId, Vector3 linear, Vector3 angular = default)
        {
            var center = GetSoftBodyCenter(softBodyId);
            TouchHost();
            foreach (int i in GetSoftBody(softBodyId).Members)
            {
                var v = linear + Vector3.Cross(angular, new Vector3(_posX[i], _posY[i], _posZ[i]) - center);
                _velX[i] = v.X;
                _velY[i] = v.Y;
                _velZ[i] = v.Z;
            }
        }

        /// <summary>
        /// Angular velocity (rad/s) of a soft body about its centre: the spin that
        /// carries the same angular momentum as its particles.
        /// </summary>
        public Vector3 GetSoftBodyAngularVelocity(int softBodyId)
        {
            SyncHost();
            var body = GetSoftBody(softBodyId);
            var sum = new float[16];
            int reference = body.Members[0];
            foreach (int i in body.Members)
            {
                float w = (_flags[i] & FlagAlive) != 0 && _invMass[i] > 0 ? 1f / _invMass[i] : 0f;
                float x = _posX[i] - _posX[reference], y = _posY[i] - _posY[reference], z = _posZ[i] - _posZ[reference];
                float vx = _velX[i], vy = _velY[i], vz = _velZ[i];
                sum[0] += w; sum[1] += w * x; sum[2] += w * y; sum[3] += w * z;
                sum[4] += w * vx; sum[5] += w * vy; sum[6] += w * vz;
                sum[7] += w * (y * vz - z * vy); sum[8] += w * (z * vx - x * vz); sum[9] += w * (x * vy - y * vx);
                sum[10] += w * x * x; sum[11] += w * x * y; sum[12] += w * x * z;
                sum[13] += w * y * y; sum[14] += w * y * z; sum[15] += w * z * z;
            }
            SolveRigidMotion(sum, out _, out _, out var omega);
            return omega;
        }

        /// <summary>
        /// Turns a soft body back into loose particles.
        /// </summary>
        public void DissolveSoftBody(int softBodyId)
        {
            TouchHost();
            var body = GetSoftBody(softBodyId);
            foreach (int i in body.Members)
                _softOf[i] = -1;
            body.Alive = false;
            body.Members = Array.Empty<int>();
            body.Clusters.Clear();
            _aliveSoftBodies--;
            _softStructureDirty = true;
        }

        private SoftBody GetSoftBody(int softBodyId)
        {
            if (!IsSoftBodyAlive(softBodyId))
                throw new ArgumentException("No such soft body.", nameof(softBodyId));
            return _softBodies[softBodyId];
        }

        #endregion

        #region Flattened Cluster Data

        /// <summary>
        /// Copies the rotations of the flat arrays back into the clusters, so that a
        /// rebuild keeps them.
        /// </summary>
        private void StoreClusterRotations()
        {
            if (_clusterCount == 0)
                return;

            int c = 0;
            foreach (var body in _softBodies)
            {
                if (!body.Alive) continue;
                foreach (var cluster in body.Clusters)
                {
                    cluster.Rotation = new Quaternion(_clusterRotation[4 * c], _clusterRotation[4 * c + 1],
                        _clusterRotation[4 * c + 2], _clusterRotation[4 * c + 3]);
                    c++;
                }
            }
        }

        private void RebuildSoftArrays()
        {
            if (!_softStructureDirty && _softBuiltForCount == _count)
                return;

            if (!_softStructureDirty)
            {
                // Only new loose particles were added: extend the per-particle table.
                ExtendSlotStart();
                return;
            }

            StoreClusterRotations();

            int clusters = 0, slots = 0;
            foreach (var body in _softBodies)
            {
                if (!body.Alive) continue;
                foreach (var cluster in body.Clusters)
                {
                    clusters++;
                    slots += cluster.Members.Length;
                }
            }

            _clusterCount = clusters;
            _slotCount = slots;
            if (_clusterStart.Length < clusters + 1) _clusterStart = new int[clusters + 1];
            if (_clusterRotation.Length < clusters * 4) _clusterRotation = new float[clusters * 4];
            if (_clusterMembers.Length < slots)
            {
                _clusterMembers = new int[slots];
                _clusterRestX = new float[slots]; _clusterRestY = new float[slots]; _clusterRestZ = new float[slots];
                _goalX = new float[slots]; _goalY = new float[slots]; _goalZ = new float[slots];
                _slots = new int[slots];
            }

            int c = 0, s = 0;
            foreach (var body in _softBodies)
            {
                if (!body.Alive) continue;
                foreach (var cluster in body.Clusters)
                {
                    _clusterStart[c] = s;
                    _clusterRotation[4 * c] = cluster.Rotation.X;
                    _clusterRotation[4 * c + 1] = cluster.Rotation.Y;
                    _clusterRotation[4 * c + 2] = cluster.Rotation.Z;
                    _clusterRotation[4 * c + 3] = cluster.Rotation.W;
                    Array.Copy(cluster.Members, 0, _clusterMembers, s, cluster.Members.Length);
                    Array.Copy(cluster.RestX, 0, _clusterRestX, s, cluster.Members.Length);
                    Array.Copy(cluster.RestY, 0, _clusterRestY, s, cluster.Members.Length);
                    Array.Copy(cluster.RestZ, 0, _clusterRestZ, s, cluster.Members.Length);
                    s += cluster.Members.Length;
                    c++;
                }
            }
            _clusterStart[clusters] = s;

            BuildDampingChunks();

            // Per particle: the slots (cluster memberships) it takes its goals from.
            ExtendSlotStart();
            _softStructureDirty = false;
            _softUploaded = false;
        }

        private void BuildDampingChunks()
        {
            int bodies = _softBodies.Count, members = 0, chunks = 0;
            foreach (var body in _softBodies)
            {
                if (!body.Alive) continue;
                members += body.Members.Length;
                chunks += (body.Members.Length + DampChunk - 1) / DampChunk;
            }

            if (_bodyStart.Length < bodies + 1) _bodyStart = new int[bodies + 1];
            if (_bodyChunkStart.Length < bodies + 1) _bodyChunkStart = new int[bodies + 1];
            if (_bodyMembers.Length < members) _bodyMembers = new int[members];
            if (_chunkBody.Length < chunks)
            {
                _chunkBody = new int[chunks];
                _chunkFrom = new int[chunks];
                _chunkTo = new int[chunks];
                _chunkPartial = new float[chunks * 16];
            }
            if (_bodyRigid.Length < bodies * 9) _bodyRigid = new float[Math.Max(bodies, 1) * 9];

            int m = 0, c = 0;
            for (int b = 0; b < bodies; b++)
            {
                _bodyStart[b] = m;
                _bodyChunkStart[b] = c;
                var body = _softBodies[b];
                if (!body.Alive) continue;

                Array.Copy(body.Members, 0, _bodyMembers, m, body.Members.Length);
                for (int from = 0; from < body.Members.Length; from += DampChunk)
                {
                    _chunkBody[c] = b;
                    _chunkFrom[c] = m + from;
                    _chunkTo[c] = m + Math.Min(from + DampChunk, body.Members.Length);
                    c++;
                }
                m += body.Members.Length;
            }
            _bodyStart[bodies] = m;
            _bodyChunkStart[bodies] = c;
            _chunkCount = c;
        }

        private void ExtendSlotStart()
        {
            int n = _count;
            if (_slotStart.Length < n + 1)
                _slotStart = new int[Math.Max(n + 1, _capacity + 1)];

            Array.Clear(_slotStart, 0, n + 1);
            for (int s = 0; s < _slotCount; s++)
                _slotStart[_clusterMembers[s] + 1]++;
            for (int i = 0; i < n; i++)
                _slotStart[i + 1] += _slotStart[i];

            var cursor = new int[n];
            Array.Copy(_slotStart, cursor, n);
            for (int s = 0; s < _slotCount; s++)
                _slots[cursor[_clusterMembers[s]]++] = s;

            _softBuiltForCount = n;
            _softUploaded = false;
        }

        /// <summary>
        /// Per-body fraction of the way to the goal covered in one sub-step. Stiffness
        /// is defined per 1/240 s, so the result does not depend on the sub-step count.
        /// </summary>
        private void ComputeSoftAlpha(float dt)
        {
            if (_softAlpha.Length < _softBodies.Count)
            {
                _softAlpha = new float[Math.Max(_softBodies.Count, 8)];
                _softBeta = new float[_softAlpha.Length];
            }

            for (int b = 0; b < _softBodies.Count; b++)
            {
                float k = _softBodies[b].Stiffness;
                _softAlpha[b] = k >= 1f ? 1f : 1f - MathF.Pow(1f - k, dt * 240f);
                float d = _softBodies[b].Damping;
                _softBeta[b] = d >= 1f ? 1f : 1f - MathF.Pow(1f - d, dt * 240f);
            }
        }

        #endregion

        #region CPU Solver

        private void MatchSoftBodies(float dt)
        {
            if (_aliveSoftBodies == 0)
                return;

            RebuildSoftArrays();
            ComputeSoftAlpha(dt);

            ForRanges(_clusterCount, (start, end) =>
            {
                for (int c = start; c < end; c++)
                    FitCluster(c);
            });

            ForRanges(_count, (start, end) =>
            {
                for (int i = start; i < end; i++)
                {
                    int from = _slotStart[i], to = _slotStart[i + 1];
                    if (to <= from)
                        continue;

                    float sx = 0, sy = 0, sz = 0;
                    for (int k = from; k < to; k++)
                    {
                        int slot = _slots[k];
                        sx += _goalX[slot]; sy += _goalY[slot]; sz += _goalZ[slot];
                    }

                    float inv = 1f / (to - from);
                    float a = _softAlpha[_softOf[i]];
                    _posX[i] += a * (sx * inv - _posX[i]);
                    _posY[i] += a * (sy * inv - _posY[i]);
                    _posZ[i] += a * (sz * inv - _posZ[i]);
                }
            });
        }

        /// <summary>
        /// Internal damping (Mueller et al., "Position Based Dynamics", 2007): the
        /// velocity of every particle is pulled towards the rigid motion of its body
        /// (the translation and rotation that carry the same momentum). Falling,
        /// sliding and spinning are untouched; squashing, stretching, wobbling and the
        /// rebound they cause are what is damped.
        /// </summary>
        private void DampSoftBodies()
        {
            if (_aliveSoftBodies == 0 || _chunkCount == 0)
                return;

            ForRanges(_chunkCount, (start, end) =>
            {
                for (int c = start; c < end; c++)
                    SumChunk(c);
            });

            for (int b = 0; b < _softBodies.Count; b++)
            {
                if (_softBodies[b].Alive && _softBeta[b] > 0)
                    SolveRigidMotion(b);
            }

            ForRanges(_count, (start, end) =>
            {
                for (int i = start; i < end; i++)
                {
                    int b = _softOf[i];
                    if (b < 0 || (_flags[i] & (FlagAlive | FlagFrozen | FlagSleeping)) != FlagAlive)
                        continue;

                    float beta = _softBeta[b];
                    if (beta <= 0)
                        continue;

                    var r = _bodyRigid;
                    float rx = _posX[i] - r[9 * b + 3], ry = _posY[i] - r[9 * b + 4], rz = _posZ[i] - r[9 * b + 5];
                    float wx = r[9 * b + 6], wy = r[9 * b + 7], wz = r[9 * b + 8];
                    float tx = r[9 * b] + (wy * rz - wz * ry);
                    float ty = r[9 * b + 1] + (wz * rx - wx * rz);
                    float tz = r[9 * b + 2] + (wx * ry - wy * rx);
                    _velX[i] += beta * (tx - _velX[i]);
                    _velY[i] += beta * (ty - _velY[i]);
                    _velZ[i] += beta * (tz - _velZ[i]);
                }
            });
        }

        /// <summary>
        /// Mass, momentum, angular momentum and second moments of one chunk of a body,
        /// relative to the body's first particle (keeps the float sums precise).
        /// </summary>
        private void SumChunk(int c)
        {
            int b = _chunkBody[c];
            int reference = _bodyMembers[_bodyStart[b]];
            float ox = _posX[reference], oy = _posY[reference], oz = _posZ[reference];

            float m = 0, mx = 0, my = 0, mz = 0, px = 0, py = 0, pz = 0, lx = 0, ly = 0, lz = 0;
            float sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0, szz = 0;
            for (int k = _chunkFrom[c]; k < _chunkTo[c]; k++)
            {
                int i = _bodyMembers[k];
                float w = (_flags[i] & FlagAlive) != 0 && _invMass[i] > 0 ? 1f / _invMass[i] : 0f;
                float x = _posX[i] - ox, y = _posY[i] - oy, z = _posZ[i] - oz;
                float vx = _velX[i], vy = _velY[i], vz = _velZ[i];
                m += w;
                mx += w * x; my += w * y; mz += w * z;
                px += w * vx; py += w * vy; pz += w * vz;
                lx += w * (y * vz - z * vy); ly += w * (z * vx - x * vz); lz += w * (x * vy - y * vx);
                sxx += w * x * x; sxy += w * x * y; sxz += w * x * z;
                syy += w * y * y; syz += w * y * z; szz += w * z * z;
            }

            var p = _chunkPartial;
            int o = 16 * c;
            p[o] = m; p[o + 1] = mx; p[o + 2] = my; p[o + 3] = mz;
            p[o + 4] = px; p[o + 5] = py; p[o + 6] = pz;
            p[o + 7] = lx; p[o + 8] = ly; p[o + 9] = lz;
            p[o + 10] = sxx; p[o + 11] = sxy; p[o + 12] = sxz; p[o + 13] = syy; p[o + 14] = syz; p[o + 15] = szz;
        }

        /// <summary>
        /// Rigid motion of a body from its chunk sums: centre-of-mass velocity, centre
        /// and angular velocity (inertia tensor inverted directly).
        /// </summary>
        private void SolveRigidMotion(int b)
        {
            var sum = new float[16];
            for (int c = _bodyChunkStart[b]; c < _bodyChunkStart[b + 1]; c++)
                for (int k = 0; k < 16; k++)
                    sum[k] += _chunkPartial[16 * c + k];

            SolveRigidMotion(sum, out var vcm, out var d, out var omega);

            int reference = _bodyMembers[_bodyStart[b]];
            var r = _bodyRigid;
            r[9 * b] = vcm.X; r[9 * b + 1] = vcm.Y; r[9 * b + 2] = vcm.Z;
            r[9 * b + 3] = _posX[reference] + d.X; r[9 * b + 4] = _posY[reference] + d.Y; r[9 * b + 5] = _posZ[reference] + d.Z;
            r[9 * b + 6] = omega.X; r[9 * b + 7] = omega.Y; r[9 * b + 8] = omega.Z;
        }

        private static void SolveRigidMotion(float[] s, out Vector3 vcm, out Vector3 d, out Vector3 omega)
        {
            float m = s[0];
            vcm = d = omega = Vector3.Zero;
            if (m <= 0)
                return;

            d = new Vector3(s[1], s[2], s[3]) / m;            // centre, relative to the reference
            vcm = new Vector3(s[4], s[5], s[6]) / m;

            // Angular momentum and second moments about the centre.
            var l = new Vector3(s[7], s[8], s[9]) - m * Vector3.Cross(d, vcm);
            float sxx = s[10] - m * d.X * d.X, sxy = s[11] - m * d.X * d.Y, sxz = s[12] - m * d.X * d.Z;
            float syy = s[13] - m * d.Y * d.Y, syz = s[14] - m * d.Y * d.Z, szz = s[15] - m * d.Z * d.Z;

            // Inertia tensor I = tr(S) E - S, slightly regularised for flat bodies.
            float trace = sxx + syy + szz;
            float eps = 1e-6f * trace + 1e-12f;
            float ixx = syy + szz + eps, iyy = sxx + szz + eps, izz = sxx + syy + eps;
            float ixy = -sxy, ixz = -sxz, iyz = -syz;

            float c00 = iyy * izz - iyz * iyz, c01 = ixz * iyz - ixy * izz, c02 = ixy * iyz - ixz * iyy;
            float det = ixx * c00 + ixy * c01 + ixz * c02;
            if (MathF.Abs(det) < 1e-30f)
                return;

            float c11 = ixx * izz - ixz * ixz, c12 = ixy * ixz - ixx * iyz, c22 = ixx * iyy - ixy * ixy;
            float inv = 1f / det;
            omega = new Vector3(
                (c00 * l.X + c01 * l.Y + c02 * l.Z) * inv,
                (c01 * l.X + c11 * l.Y + c12 * l.Z) * inv,
                (c02 * l.X + c12 * l.Y + c22 * l.Z) * inv);
        }

        private void FitCluster(int c)
        {
            int from = _clusterStart[c], to = _clusterStart[c + 1];

            float total = 0, cx = 0, cy = 0, cz = 0;
            for (int k = from; k < to; k++)
            {
                int i = _clusterMembers[k];
                float m = (_flags[i] & FlagAlive) != 0 && _invMass[i] > 0 ? 1f / _invMass[i] : 0f;
                cx += m * _posX[i]; cy += m * _posY[i]; cz += m * _posZ[i];
                total += m;
            }
            if (total <= 0)
                return;
            cx /= total; cy /= total; cz /= total;

            Vector3 a0 = default, a1 = default, a2 = default;
            for (int k = from; k < to; k++)
            {
                int i = _clusterMembers[k];
                float m = (_flags[i] & FlagAlive) != 0 && _invMass[i] > 0 ? 1f / _invMass[i] : 0f;
                var d = m * new Vector3(_posX[i] - cx, _posY[i] - cy, _posZ[i] - cz);
                a0 += d * _clusterRestX[k];
                a1 += d * _clusterRestY[k];
                a2 += d * _clusterRestZ[k];
            }

            var q = new Quaternion(_clusterRotation[4 * c], _clusterRotation[4 * c + 1], _clusterRotation[4 * c + 2], _clusterRotation[4 * c + 3]);
            q = ExtractRotation(a0, a1, a2, q);
            _clusterRotation[4 * c] = q.X;
            _clusterRotation[4 * c + 1] = q.Y;
            _clusterRotation[4 * c + 2] = q.Z;
            _clusterRotation[4 * c + 3] = q.W;

            var center = new Vector3(cx, cy, cz);
            for (int k = from; k < to; k++)
            {
                var g = center + Vector3.Transform(new Vector3(_clusterRestX[k], _clusterRestY[k], _clusterRestZ[k]), q);
                _goalX[k] = g.X;
                _goalY[k] = g.Y;
                _goalZ[k] = g.Z;
            }
        }

        #endregion

        #region GPU Glue

#if !NETSTANDARD2_1
        private void UploadSoftBodies()
        {
            if (_aliveSoftBodies == 0)
                return;

            RebuildSoftArrays();
            if (_softUploaded)
                return;

            var g = _gpuCompute!;
            var rt = g.Runtime!;
            g.EnsureSoftCapacity(Math.Max(_clusterCount, 1), Math.Max(_slotCount, 1), Math.Max(_softBodies.Count, 1));

            rt.Write(g.ClusterStart, _clusterStart, _clusterCount + 1);
            rt.Write(g.Members, _clusterMembers, _slotCount);
            rt.Write(g.RestX, _clusterRestX, _slotCount);
            rt.Write(g.RestY, _clusterRestY, _slotCount);
            rt.Write(g.RestZ, _clusterRestZ, _slotCount);
            rt.Write(g.Rotation, _clusterRotation, _clusterCount * 4);
            rt.Write(g.SlotStart, _slotStart, _count + 1);
            rt.Write(g.Slots, _slots, _slotCount);
            if (_chunkCount > 0)
            {
                g.EnsureDampCapacity(_chunkCount, _bodyStart[_softBodies.Count], _softBodies.Count);
                rt.Write(g.BodyStart, _bodyStart, _softBodies.Count + 1);
                rt.Write(g.BodyMembers, _bodyMembers, _bodyStart[_softBodies.Count]);
                rt.Write(g.BodyChunkStart, _bodyChunkStart, _softBodies.Count + 1);
                rt.Write(g.ChunkBody, _chunkBody, _chunkCount);
                rt.Write(g.ChunkFrom, _chunkFrom, _chunkCount);
                rt.Write(g.ChunkTo, _chunkTo, _chunkCount);
            }
            _softUploaded = true;
        }

        private void StepSoftBodiesGpu(MassiveGpuCompute g, OpenCLRuntime rt, int n, float dt)
        {
            if (_aliveSoftBodies == 0 || _clusterCount == 0)
                return;

            ComputeSoftAlpha(dt);
            rt.Write(g.Alpha, _softAlpha, _softBodies.Count);
            rt.Write(g.Beta, _softBeta, _softBodies.Count);

            rt.Run("soft_fit", _clusterCount, LocalSize, _clusterCount, g.ClusterStart, g.Members,
                g.RestX, g.RestY, g.RestZ, g.Rotation, g.PX, g.PY, g.PZ, g.InvMass, g.Flags,
                g.GoalX, g.GoalY, g.GoalZ);
            rt.Run("soft_apply", n, LocalSize, n, g.SlotStart, g.Slots, g.Body, g.Alpha,
                g.GoalX, g.GoalY, g.GoalZ, g.PX, g.PY, g.PZ);
        }

        private void DampSoftBodiesGpu(MassiveGpuCompute g, OpenCLRuntime rt, int n)
        {
            if (_aliveSoftBodies == 0 || _chunkCount == 0)
                return;

            int bodies = _softBodies.Count;
            rt.Run("soft_body_sum", _chunkCount, LocalSize, _chunkCount, g.ChunkBody, g.ChunkFrom, g.ChunkTo,
                g.BodyStart, g.BodyMembers, g.PX, g.PY, g.PZ, g.VX, g.VY, g.VZ, g.InvMass, g.Flags, g.ChunkPartial);
            rt.Run("soft_body_rigid", bodies, LocalSize, bodies, g.BodyChunkStart, g.ChunkPartial,
                g.BodyStart, g.BodyMembers, g.PX, g.PY, g.PZ, g.Beta, g.BodyRigid);
            rt.Run("soft_body_damp", n, LocalSize, n, g.Body, g.Beta, g.BodyRigid, g.Flags,
                g.PX, g.PY, g.PZ, g.VX, g.VY, g.VZ);
        }

        private void DownloadSoftBodies()
        {
            if (_aliveSoftBodies == 0 || _clusterCount == 0 || !_softUploaded)
                return;

            var g = _gpuCompute!;
            g.Runtime!.Read(g.Rotation, _clusterRotation, _clusterCount * 4);
        }
#endif

        #endregion

        #region Maintenance

        private void ClearSoftBodies()
        {
            _softBodies.Clear();
            _aliveSoftBodies = 0;
            for (int i = 0; i < _softOf.Length; i++)
                _softOf[i] = -1;
            _clusterCount = 0;
            _slotCount = 0;
            _softStructureDirty = true;
        }

        /// <summary>
        /// Rewrites the soft bodies after <see cref="Compact"/> has moved particles.
        /// Removed particles leave their clusters; the remaining rest shape is
        /// re-centred.
        /// </summary>
        private void RemapSoftBodies(int[] newIndex)
        {
            StoreClusterRotations();

            foreach (var body in _softBodies)
            {
                if (!body.Alive)
                    continue;

                var kept = new List<int>(body.Members.Length);
                foreach (int i in body.Members)
                    if (newIndex[i] >= 0) kept.Add(newIndex[i]);
                body.Members = kept.ToArray();

                var clusters = new List<Cluster>();
                foreach (var cluster in body.Clusters)
                {
                    var members = new List<int>();
                    var rx = new List<float>(); var ry = new List<float>(); var rz = new List<float>();
                    for (int k = 0; k < cluster.Members.Length; k++)
                    {
                        int moved = newIndex[cluster.Members[k]];
                        if (moved < 0) continue;
                        members.Add(moved);
                        rx.Add(cluster.RestX[k]); ry.Add(cluster.RestY[k]); rz.Add(cluster.RestZ[k]);
                    }

                    if (members.Count < 4)
                        continue;

                    // Re-centre the rest shape on the members that are left.
                    float mx = 0, my = 0, mz = 0;
                    for (int k = 0; k < members.Count; k++) { mx += rx[k]; my += ry[k]; mz += rz[k]; }
                    mx /= members.Count; my /= members.Count; mz /= members.Count;
                    for (int k = 0; k < members.Count; k++) { rx[k] -= mx; ry[k] -= my; rz[k] -= mz; }

                    clusters.Add(new Cluster
                    {
                        Members = members.ToArray(),
                        RestX = rx.ToArray(), RestY = ry.ToArray(), RestZ = rz.ToArray(),
                        Rotation = cluster.Rotation
                    });
                }

                body.Clusters = clusters;
                if (body.Members.Length == 0 || clusters.Count == 0)
                {
                    foreach (int i in body.Members) _softOf[i] = -1;
                    body.Alive = false;
                    _aliveSoftBodies--;
                }
            }

            _softStructureDirty = true;
        }

        #endregion
    }
}
