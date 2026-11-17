using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;
using Quaternion = System.Numerics.Quaternion;

namespace Artemis.Particles
{
    // Meshes made of particles.
    //
    // A mesh can be turned into particles in two ways:
    //   * frozen grains that hold the shape until they are released and then behave
    //     like sand (SpawnFromMesh), or
    //   * a rigid group: the particles move together as one solid object, which
    //     shatters into smaller rigid fragments and loose grains when it is hit
    //     hard enough (SpawnRigidBodyFromMesh) - a glass falling on the floor.
    //
    // Rigid groups use shape matching: every sub-step the particles are first moved
    // individually (gravity, contacts), then the rigid placement that fits their new
    // positions best is computed and the particles are snapped onto it.
    public partial class MassiveParticleSystem
    {
        #region Types and Fields

        private sealed class RigidGroup
        {
            public int[] Members = Array.Empty<int>();
            public Vector3 Center;
            public Quaternion Rotation = Quaternion.Identity;
            public Vector3 Velocity;
            public float BreakSpeed;
            public int Generation;
            public bool Alive = true;

            // Result of the last fit, consumed by the shatter check.
            public float ImpactSpeed;
            public Vector3 ImpactPoint;
        }

        private readonly List<RigidGroup> _groups = new();
        private int[] _groupOf = Array.Empty<int>();
        private float[] _restX = Array.Empty<float>();
        private float[] _restY = Array.Empty<float>();
        private float[] _restZ = Array.Empty<float>();
        private int _aliveGroups;
        private Random _shatterRandom = new(12345);

        #endregion

        #region Properties and Events

        /// <summary>
        /// Number of rigid groups that currently exist (fragments included).
        /// </summary>
        public int RigidGroupCount => _aliveGroups;

        /// <summary>
        /// Number of fragments a rigid group breaks into when it shatters.
        /// </summary>
        public int ShatterFragments { get; set; } = 24;

        /// <summary>
        /// Fragments with fewer particles than this are not kept as rigid pieces but
        /// crumble into loose grains.
        /// </summary>
        public int MinFragmentParticles { get; set; } = 6;

        /// <summary>
        /// Factor by which fragments are tougher than the object they came from, so
        /// that one impact does not grind everything to dust at once.
        /// </summary>
        public float FragmentToughness { get; set; } = 2.5f;

        /// <summary>
        /// How many times a piece may shatter again. 0 means the original object
        /// breaks once and its fragments are unbreakable.
        /// </summary>
        public int MaxShatterDepth { get; set; } = 2;

        /// <summary>
        /// Speed (m/s) with which fragments are pushed away from the point of impact.
        /// </summary>
        public float ShatterScatterSpeed { get; set; } = 1.0f;

        /// <summary>
        /// Raised when a rigid group shatters: group id, point of impact and the
        /// number of rigid fragments created.
        /// </summary>
        public event Action<int, Vector3, int>? GroupShattered;

        #endregion

        #region Allocation

        private void AllocateGroupArrays(int capacity)
        {
            _groupOf = new int[capacity];
            for (int i = 0; i < capacity; i++)
                _groupOf[i] = -1;

            _restX = new float[capacity];
            _restY = new float[capacity];
            _restZ = new float[capacity];
        }

        #endregion

        #region Spawning from Meshes

        /// <summary>
        /// Rebuilds a mesh out of particles.
        /// </summary>
        /// <param name="mesh">The mesh, already placed where it should appear.</param>
        /// <param name="particleRadius">Radius of the grains; they are spaced one diameter apart.</param>
        /// <param name="firstIndex">Index of the first particle created.</param>
        /// <param name="solid">Fill the volume (closed meshes) or only cover the surface.</param>
        /// <param name="frozen">
        /// Keep the grains in place until they are released with <see cref="ReleaseAll"/>,
        /// <see cref="ReleaseInRadius"/>, <see cref="ReleaseBeyondPlane"/> or an impact
        /// above <see cref="FrozenBreakSpeed"/>.
        /// </param>
        /// <param name="mass">Mass of one grain.</param>
        /// <param name="color">Colour of the grains.</param>
        /// <param name="jitter">
        /// Random displacement of each grain as a fraction of the spacing. A perfect
        /// lattice of spheres can stay standing like a house of cards after it is
        /// released; a little disorder makes it slump the way sand does.
        /// </param>
        /// <param name="seed">Seed of the jitter.</param>
        /// <returns>The number of particles created (limited by the remaining capacity).</returns>
        public int SpawnFromMesh(
            TriangleMesh mesh,
            float particleRadius,
            out int firstIndex,
            bool solid = true,
            bool frozen = true,
            float mass = 1.0f,
            uint color = 0xFFFFFFFF,
            float jitter = 0.15f,
            int seed = 12345)
        {
            if (mesh == null) throw new ArgumentNullException(nameof(mesh));
            if (!(particleRadius > 0)) throw new ArgumentOutOfRangeException(nameof(particleRadius));

            firstIndex = _count;
            int room = _capacity - _count;
            if (room <= 0)
                return 0;

            float spacing = 2f * particleRadius;
            var points = solid
                ? MeshSampler.FillVolume(mesh, spacing, room)
                : MeshSampler.SampleSurface(mesh, spacing, room);

            if (points.Count == 0)
                return 0;

            if (jitter > 0)
            {
                var random = new Random(seed);
                float amplitude = jitter * spacing;
                for (int i = 0; i < points.Count; i++)
                {
                    points[i] += new Vector3(
                        ((float)random.NextDouble() - 0.5f) * amplitude,
                        ((float)random.NextDouble() - 0.5f) * amplitude,
                        ((float)random.NextDouble() - 0.5f) * amplitude);
                }
            }

            int created = SpawnBatch(points, mass, particleRadius, float.MaxValue, color);
            if (created > 0 && frozen)
                Freeze(firstIndex, created);

            return created;
        }

        /// <summary>
        /// Rebuilds a mesh out of particles that move together as one rigid object.
        /// The object shatters into fragments when an impact changes its velocity by
        /// more than <paramref name="breakSpeed"/>.
        /// </summary>
        /// <param name="mesh">The mesh, already placed where it should appear.</param>
        /// <param name="particleRadius">
        /// Radius of the particles. For thin-walled shapes (a glass, a vase) it must
        /// not exceed half the wall thickness, or the walls come out with holes.
        /// </param>
        /// <param name="breakSpeed">
        /// Change of velocity (m/s) in one sub-step that breaks the object. A drop
        /// from height h lands at sqrt(2 g h): about 4.4 m/s from one metre.
        /// </param>
        /// <param name="groupId">Id of the rigid group, or -1 if nothing was created.</param>
        /// <param name="solid">Fill the volume or only cover the surface.</param>
        /// <param name="mass">Mass of one particle.</param>
        /// <param name="color">Colour of the particles.</param>
        /// <returns>The number of particles created.</returns>
        public int SpawnRigidBodyFromMesh(
            TriangleMesh mesh,
            float particleRadius,
            float breakSpeed,
            out int groupId,
            bool solid = true,
            float mass = 1.0f,
            uint color = 0xFFFFFFFF)
        {
            int created = SpawnFromMesh(mesh, particleRadius, out int first, solid, frozen: false, mass, color, jitter: 0);
            groupId = created > 0 ? CreateRigidGroup(first, created, breakSpeed) : -1;
            return created;
        }

        #endregion

        #region Rigid Groups

        /// <summary>
        /// Binds a range of particles into one rigid object.
        /// </summary>
        /// <param name="firstIndex">First particle of the range.</param>
        /// <param name="count">Number of particles.</param>
        /// <param name="breakSpeed">
        /// Change of velocity (m/s) in one sub-step that shatters the object;
        /// infinity makes it unbreakable.
        /// </param>
        /// <returns>The id of the group.</returns>
        public int CreateRigidGroup(int firstIndex, int count, float breakSpeed = float.PositiveInfinity)
        {
            if (firstIndex < 0 || count <= 0 || firstIndex + count > _count)
                throw new ArgumentOutOfRangeException(nameof(count));

            var members = new int[count];
            for (int i = 0; i < count; i++)
                members[i] = firstIndex + i;

            return AddGroup(members, breakSpeed, generation: 0);
        }

        private int AddGroup(int[] members, float breakSpeed, int generation)
        {
            var group = new RigidGroup
            {
                Members = members,
                BreakSpeed = breakSpeed,
                Generation = generation
            };

            int id = _groups.Count;
            _groups.Add(group);
            _aliveGroups++;

            // Centre of mass and rest shape, taken from where the particles are now.
            double cx = 0, cy = 0, cz = 0, vx = 0, vy = 0, vz = 0, total = 0;
            foreach (int i in members)
            {
                double m = _mass[i];
                cx += m * _posX[i]; cy += m * _posY[i]; cz += m * _posZ[i];
                vx += m * _velX[i]; vy += m * _velY[i]; vz += m * _velZ[i];
                total += m;
            }

            if (total <= 0) total = 1;
            group.Center = new Vector3((float)(cx / total), (float)(cy / total), (float)(cz / total));
            group.Velocity = new Vector3((float)(vx / total), (float)(vy / total), (float)(vz / total));

            foreach (int i in members)
            {
                _groupOf[i] = id;
                _restX[i] = _posX[i] - group.Center.X;
                _restY[i] = _posY[i] - group.Center.Y;
                _restZ[i] = _posZ[i] - group.Center.Z;

                // Members of a rigid object are neither frozen nor put to sleep one by one.
                _flags[i] &= unchecked((byte)~(FlagFrozen | FlagSleeping));
                _rest[i] = 0;
            }

            return id;
        }

        /// <summary>
        /// Position and orientation of a rigid group.
        /// </summary>
        public (Vector3 Center, Quaternion Rotation) GetGroupPose(int groupId)
        {
            var group = GetGroup(groupId);
            return (group.Center, group.Rotation);
        }

        /// <summary>
        /// Whether a rigid group still exists (it has not shattered or been dissolved).
        /// </summary>
        public bool IsGroupAlive(int groupId)
            => groupId >= 0 && groupId < _groups.Count && _groups[groupId].Alive;

        /// <summary>
        /// Number of particles in a rigid group.
        /// </summary>
        public int GetGroupParticleCount(int groupId) => GetGroup(groupId).Members.Length;

        /// <summary>
        /// Id of the rigid group a particle belongs to, or -1 for a loose particle.
        /// </summary>
        public int GetParticleGroup(int index) => _groupOf[index];

        /// <summary>
        /// Sets the linear and angular velocity of a rigid group (throw it, spin it).
        /// </summary>
        public void SetGroupVelocity(int groupId, Vector3 linear, Vector3 angular = default)
        {
            var group = GetGroup(groupId);
            group.Velocity = linear;

            foreach (int i in group.Members)
            {
                var offset = new Vector3(_posX[i], _posY[i], _posZ[i]) - group.Center;
                var v = linear + Vector3.Cross(angular, offset);
                _velX[i] = v.X;
                _velY[i] = v.Y;
                _velZ[i] = v.Z;
            }
        }

        /// <summary>
        /// Turns a rigid group back into loose particles.
        /// </summary>
        public void DissolveGroup(int groupId)
        {
            var group = GetGroup(groupId);
            foreach (int i in group.Members)
                _groupOf[i] = -1;

            group.Alive = false;
            group.Members = Array.Empty<int>();
            _aliveGroups--;
        }

        /// <summary>
        /// Shatters a rigid group around a point, regardless of its break speed.
        /// </summary>
        /// <returns>The number of rigid fragments created.</returns>
        public int ShatterGroup(int groupId, Vector3 impactPoint)
            => Shatter(groupId, impactPoint, 0);

        private RigidGroup GetGroup(int groupId)
        {
            if (groupId < 0 || groupId >= _groups.Count || !_groups[groupId].Alive)
                throw new ArgumentException("No such rigid group.", nameof(groupId));
            return _groups[groupId];
        }

        #endregion

        #region Shape Matching

        /// <summary>
        /// Snaps the particles of every rigid group onto the rigid placement that best
        /// fits where the sub-step has moved them, and shatters the groups that were
        /// stopped too abruptly.
        /// </summary>
        private void MatchRigidGroups(float dt)
        {
            if (_aliveGroups == 0)
                return;

            int groupCount = _groups.Count;
            if (groupCount >= 8 && _threadCount > 1)
                Parallel.For(0, groupCount, _parallelOptions, g => FitGroup(_groups[g], dt));
            else
                for (int g = 0; g < groupCount; g++) FitGroup(_groups[g], dt);

            // Shattering changes the group list, so it happens after the parallel fit.
            for (int g = 0; g < groupCount; g++)
            {
                var group = _groups[g];
                if (group.Alive && group.ImpactSpeed > group.BreakSpeed)
                    Shatter(g, group.ImpactPoint, dt);
            }
        }

        private void FitGroup(RigidGroup group, float dt)
        {
            group.ImpactSpeed = 0;
            if (!group.Alive)
                return;

            var members = group.Members;
            float minX = (float)Bounds.Min.X, minY = (float)Bounds.Min.Y, minZ = (float)Bounds.Min.Z;
            float maxX = (float)Bounds.Max.X, maxY = (float)Bounds.Max.Y, maxZ = (float)Bounds.Max.Z;

            // Centre of mass of the moved particles. Particles pushed into a wall are
            // clamped first: that is what makes the object feel the wall and tip over.
            double cx = 0, cy = 0, cz = 0, px = 0, py = 0, pz = 0, total = 0;
            double hitX = 0, hitY = 0, hitZ = 0, hitWeight = 0;

            foreach (int i in members)
            {
                float r = _radius[i];
                float x = Math.Clamp(_posX[i], minX + r, Math.Max(minX + r, maxX - r));
                float y = Math.Clamp(_posY[i], minY + r, Math.Max(minY + r, maxY - r));
                float z = Math.Clamp(_posZ[i], minZ + r, Math.Max(minZ + r, maxZ - r));

                // Particles that were pushed back (by a wall here, or by other
                // particles in the contact pass) mark where the object was hit.
                float pushed = Math.Abs(x - _posX[i]) + Math.Abs(y - _posY[i]) + Math.Abs(z - _posZ[i])
                               + Math.Abs(_dpX[i]) + Math.Abs(_dpY[i]) + Math.Abs(_dpZ[i]);
                if (pushed > 0)
                {
                    hitX += pushed * x; hitY += pushed * y; hitZ += pushed * z;
                    hitWeight += pushed;
                }

                _posX[i] = x; _posY[i] = y; _posZ[i] = z;

                double m = _mass[i];
                cx += m * x; cy += m * y; cz += m * z;
                px += m * _prevX[i]; py += m * _prevY[i]; pz += m * _prevZ[i];
                total += m;
            }

            if (total <= 0)
                return;

            var center = new Vector3((float)(cx / total), (float)(cy / total), (float)(cz / total));
            var previousCenter = new Vector3((float)(px / total), (float)(py / total), (float)(pz / total));

            // Covariance between the current offsets and the rest shape.
            double a00 = 0, a01 = 0, a02 = 0, a10 = 0, a11 = 0, a12 = 0, a20 = 0, a21 = 0, a22 = 0;
            foreach (int i in members)
            {
                double m = _mass[i];
                double dx = _posX[i] - center.X, dy = _posY[i] - center.Y, dz = _posZ[i] - center.Z;
                double qx = _restX[i], qy = _restY[i], qz = _restZ[i];
                a00 += m * dx * qx; a01 += m * dx * qy; a02 += m * dx * qz;
                a10 += m * dy * qx; a11 += m * dy * qy; a12 += m * dy * qz;
                a20 += m * dz * qx; a21 += m * dz * qy; a22 += m * dz * qz;
            }

            var column0 = new Vector3((float)a00, (float)a10, (float)a20);
            var column1 = new Vector3((float)a01, (float)a11, (float)a21);
            var column2 = new Vector3((float)a02, (float)a12, (float)a22);
            var rotation = ExtractRotation(column0, column1, column2, group.Rotation);

            // Particles are moved along straight lines, so a body turning by an angle
            // a per sub-step is fitted with atan(a) instead: slightly too little, and
            // the spin would slowly die. Undo that for the rotation of this sub-step.
            var delta = Quaternion.Normalize(rotation * Quaternion.Conjugate(group.Rotation));
            if (delta.W < 0) delta = Quaternion.Negate(delta);
            float sinHalf = MathF.Sqrt(delta.X * delta.X + delta.Y * delta.Y + delta.Z * delta.Z);
            if (sinHalf > 1e-7f)
            {
                float angle = 2f * MathF.Atan2(sinHalf, delta.W);
                if (angle < 1.0f)
                {
                    var axis = new Vector3(delta.X, delta.Y, delta.Z) / sinHalf;
                    rotation = Quaternion.Normalize(
                        Quaternion.CreateFromAxisAngle(axis, MathF.Tan(angle)) * group.Rotation);
                }
            }

            // Goal positions.
            float lowX = 0, lowY = 0, lowZ = 0, highX = 0, highY = 0, highZ = 0;
            foreach (int i in members)
            {
                var goal = center + Vector3.Transform(new Vector3(_restX[i], _restY[i], _restZ[i]), rotation);
                _posX[i] = goal.X; _posY[i] = goal.Y; _posZ[i] = goal.Z;

                float r = _radius[i];
                lowX = Math.Max(lowX, minX + r - goal.X); highX = Math.Max(highX, goal.X - (maxX - r));
                lowY = Math.Max(lowY, minY + r - goal.Y); highY = Math.Max(highY, goal.Y - (maxY - r));
                lowZ = Math.Max(lowZ, minZ + r - goal.Z); highZ = Math.Max(highZ, goal.Z - (maxZ - r));
            }

            // The rigid fit can leave part of the object inside a wall. Move the whole
            // object out; otherwise it would slowly sink, because only the particles
            // touching the wall would be held back.
            var shift = new Vector3(lowX - highX, lowY - highY, lowZ - highZ);
            if (shift != Vector3.Zero)
            {
                // Coulomb friction against the wall on the object as a whole.
                var move = center + shift - previousCenter;
                float depth = shift.Length();
                var normal = shift / depth;
                var tangential = move - normal * Vector3.Dot(move, normal);
                float slide = tangential.Length();
                if (slide > 1e-9f && Friction > 0)
                    shift -= tangential * Math.Min(1f, Friction * depth / slide);

                foreach (int i in members)
                {
                    _posX[i] += shift.X; _posY[i] += shift.Y; _posZ[i] += shift.Z;
                }
                center += shift;
            }

            // The correction of each member is now part of the rigid motion.
            foreach (int i in members)
            {
                _stabX[i] = 0; _stabY[i] = 0; _stabZ[i] = 0;
            }

            // How abruptly was the object stopped? Compare the velocity it ends up
            // with to the one it would have had in free flight.
            var velocity = (center - previousCenter) / dt;
            var expected = (group.Velocity + Gravity * dt) * Damping;
            group.ImpactSpeed = (velocity - expected).Length();
            group.ImpactPoint = hitWeight > 0
                ? new Vector3((float)(hitX / hitWeight), (float)(hitY / hitWeight), (float)(hitZ / hitWeight))
                : center;

            group.Center = center;
            group.Rotation = rotation;
            group.Velocity = velocity;
        }

        /// <summary>
        /// Rotation closest to a 3x3 matrix given by its columns, found by iterating
        /// from the previous rotation (Mueller et al., "A robust method to extract the
        /// rotational part of deformations", 2016). Stable even for flat or
        /// degenerate particle sets, where a polar decomposition breaks down.
        /// </summary>
        private static Quaternion ExtractRotation(Vector3 a0, Vector3 a1, Vector3 a2, Quaternion q)
        {
            for (int iteration = 0; iteration < 20; iteration++)
            {
                var r0 = Vector3.Transform(Vector3.UnitX, q);
                var r1 = Vector3.Transform(Vector3.UnitY, q);
                var r2 = Vector3.Transform(Vector3.UnitZ, q);

                var omega = Vector3.Cross(r0, a0) + Vector3.Cross(r1, a1) + Vector3.Cross(r2, a2);
                float scale = Math.Abs(Vector3.Dot(r0, a0) + Vector3.Dot(r1, a1) + Vector3.Dot(r2, a2)) + 1e-12f;
                omega /= scale;

                float angle = omega.Length();
                if (angle < 1e-7f)
                    break;

                q = Quaternion.Normalize(Quaternion.CreateFromAxisAngle(omega / angle, angle) * q);
            }

            return q;
        }

        #endregion

        #region Shattering

        private int Shatter(int groupId, Vector3 impactPoint, float dt)
        {
            var group = GetGroup(groupId);
            var parentVelocity = group.Velocity;
            var members = group.Members;
            int fragmentsWanted = Math.Max(2, Math.Min(ShatterFragments, members.Length));

            // Fragment seeds are members of the object, drawn so that they crowd
            // around the point of impact: small shards there, large pieces far away.
            float size = 0;
            foreach (int i in members)
            {
                float d = (new Vector3(_posX[i], _posY[i], _posZ[i]) - impactPoint).Length();
                if (d > size) size = d;
            }
            if (size <= 0) size = 1;

            var seeds = new Vector3[fragmentsWanted];
            for (int s = 0; s < fragmentsWanted; s++)
            {
                int best = members[_shatterRandom.Next(members.Length)];
                float bestScore = float.MaxValue;

                // A few random candidates, keeping the one closest to the impact.
                for (int attempt = 0; attempt < 3; attempt++)
                {
                    int candidate = members[_shatterRandom.Next(members.Length)];
                    float d = (new Vector3(_posX[candidate], _posY[candidate], _posZ[candidate]) - impactPoint).Length();
                    float score = d / size + (float)_shatterRandom.NextDouble() * 0.5f;
                    if (score < bestScore)
                    {
                        bestScore = score;
                        best = candidate;
                    }
                }

                seeds[s] = new Vector3(_posX[best], _posY[best], _posZ[best]);
            }

            // Every particle joins the nearest seed (a Voronoi partition of the object).
            var owner = new int[members.Length];
            var sizes = new int[fragmentsWanted];
            Parallel.For(0, members.Length, _parallelOptions, k =>
            {
                int i = members[k];
                var p = new Vector3(_posX[i], _posY[i], _posZ[i]);
                int nearest = 0;
                float nearestSq = float.MaxValue;
                for (int s = 0; s < seeds.Length; s++)
                {
                    float d = Vector3.DistanceSquared(p, seeds[s]);
                    if (d < nearestSq)
                    {
                        nearestSq = d;
                        nearest = s;
                    }
                }
                owner[k] = nearest;
            });

            foreach (int o in owner)
                sizes[o]++;

            var fragments = new int[fragmentsWanted][];
            var fill = new int[fragmentsWanted];
            for (int s = 0; s < fragmentsWanted; s++)
                fragments[s] = new int[sizes[s]];
            for (int k = 0; k < members.Length; k++)
                fragments[owner[k]][fill[owner[k]]++] = members[k];

            // The original object is gone.
            int generation = group.Generation + 1;
            float fragmentBreakSpeed = generation > MaxShatterDepth
                ? float.PositiveInfinity
                : group.BreakSpeed * FragmentToughness;

            foreach (int i in members)
                _groupOf[i] = -1;
            group.Alive = false;
            group.Members = Array.Empty<int>();
            _aliveGroups--;

            int created = 0;
            for (int s = 0; s < fragmentsWanted; s++)
            {
                var piece = fragments[s];
                if (piece.Length == 0)
                    continue;

                // Push the piece away from the impact.
                var pieceCenter = Vector3.Zero;
                foreach (int i in piece)
                    pieceCenter += new Vector3(_posX[i], _posY[i], _posZ[i]);
                pieceCenter /= piece.Length;

                var away = pieceCenter - impactPoint;
                float distance = away.Length();
                var kick = distance > 1e-6f
                    ? away / distance * (ShatterScatterSpeed * (0.5f + (float)_shatterRandom.NextDouble()))
                    : Vector3.Zero;

                foreach (int i in piece)
                {
                    if (dt > 0)
                    {
                        // In the middle of a sub-step velocities are derived from the
                        // distance travelled, so the kick is given as a displacement.
                        _posX[i] += kick.X * dt;
                        _posY[i] += kick.Y * dt;
                        _posZ[i] += kick.Z * dt;
                    }
                    else
                    {
                        _velX[i] += kick.X;
                        _velY[i] += kick.Y;
                        _velZ[i] += kick.Z;
                    }
                }

                // Tiny pieces become loose grains (dust).
                if (piece.Length < MinFragmentParticles)
                    continue;

                int id = AddGroup(piece, fragmentBreakSpeed, generation);
                if (dt > 0)
                    _groups[id].Velocity = parentVelocity + kick;
                created++;
            }

            GroupShattered?.Invoke(groupId, impactPoint, created);
            return created;
        }

        #endregion

        #region Maintenance

        private void ClearGroups()
        {
            _groups.Clear();
            _aliveGroups = 0;
            for (int i = 0; i < _groupOf.Length; i++)
                _groupOf[i] = -1;
        }

        /// <summary>
        /// Rewrites the member lists after <see cref="Compact"/> has moved particles.
        /// </summary>
        /// <param name="newIndex">New index of every old particle, or -1 if it was removed.</param>
        private void RemapGroups(int[] newIndex)
        {
            foreach (var group in _groups)
            {
                if (!group.Alive)
                    continue;

                var kept = new List<int>(group.Members.Length);
                foreach (int i in group.Members)
                {
                    if (newIndex[i] >= 0)
                        kept.Add(newIndex[i]);
                }
                group.Members = kept.ToArray();

                if (group.Members.Length == 0)
                {
                    group.Alive = false;
                    _aliveGroups--;
                }
            }
        }

        #endregion
    }
}
