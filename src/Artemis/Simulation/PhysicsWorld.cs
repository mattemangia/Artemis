using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using Artemis.Bodies;
using Artemis.Collision;
using Artemis.Core;
using Artemis.Forces;
using Artemis.Particles;

namespace Artemis.Simulation
{
    /// <summary>
    /// The main physics simulation world with multi-threading support.
    /// Manages bodies, forces, and updates with parallel collision detection.
    /// </summary>
    public class PhysicsWorld
    {
        #region Fields

        private readonly List<IPhysicsBody> _bodies;
        private readonly List<IForce> _globalForces;
        private readonly List<ParticleSystem> _particleSystems;
        private readonly List<CollisionInfo> _collisions;
        private double _accumulator;

        // Spatial hashing for broad phase
        private readonly Dictionary<long, List<int>> _spatialGrid;
        private readonly Stack<List<int>> _cellListPool = new();
        private readonly List<int> _largeBodies = new();
        private readonly List<IPhysicsBody> _sweepCandidates = new();
        private readonly List<(int A, int B)> _pairs = new();
        private (int, int, int)[] _cellMin = new (int, int, int)[64];
        private (int, int, int)[] _cellMax = new (int, int, int)[64];
        private CollisionInfo[] _narrowPhaseResults = new CollisionInfo[64];
        private const double MaxCellsPerBody = 64;
        private const double CellCoordLimit = 1 << 20;
        private const int ParallelThreshold = 256;

        // Contact solving and sleeping
        private readonly ContactSolver _solver = new();
        private readonly SleepManager _sleep = new();
        private volatile IPhysicsBody[]? _bodySnapshot;
        private Vector3D _gravity = new(0, -PhysicsConstants.EarthGravity, 0);
        private double _broadPhaseCellSize = 2.0;

        // Thread synchronization
        private readonly ReaderWriterLockSlim _bodiesLock = new();

        #endregion

        #region Properties

        /// <summary>
        /// Gets all bodies in the world.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> Bodies
        {
            get => GetBodySnapshot();
        }

        /// <summary>
        /// Gets all global forces.
        /// </summary>
        public IReadOnlyList<IForce> GlobalForces => _globalForces;

        /// <summary>
        /// Gets all particle systems.
        /// </summary>
        public IReadOnlyList<ParticleSystem> ParticleSystems => _particleSystems;

        /// <summary>
        /// Gets the collisions from the last update.
        /// </summary>
        public IReadOnlyList<CollisionInfo> LastCollisions => _collisions;

        /// <summary>
        /// Gets or sets the fixed time step for physics simulation.
        /// </summary>
        public double FixedTimeStep { get; set; } = PhysicsConstants.DefaultFixedTimeStep;

        /// <summary>
        /// Gets or sets the maximum number of sub-steps per frame.
        /// </summary>
        public int MaxSubSteps { get; set; } = PhysicsConstants.DefaultMaxSubSteps;

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity
        {
            get => _gravity;
            set
            {
                if (_gravity == value)
                    return;

                _gravity = value;
                // Sleeping bodies must notice that "down" has changed.
                WakeAllBodies();
            }
        }

        /// <summary>
        /// Gets or sets the number of velocity iterations of the contact solver.
        /// Higher values make stacks and piles stiffer at a linear cost.
        /// </summary>
        public int SolverIterations { get; set; } = 12;

        /// <summary>
        /// Gets or sets the number of penetration recovery iterations.
        /// </summary>
        public int PositionIterations { get; set; } = 4;

        /// <summary>
        /// Gets or sets whether resting bodies may be put to sleep.
        /// </summary>
        public bool AllowSleeping { get; set; } = true;

        /// <summary>
        /// Linear speed (m/s) below which a body counts as resting.
        /// </summary>
        public double SleepVelocityThreshold
        {
            get => _sleep.SleepVelocityThreshold;
            set => _sleep.SleepVelocityThreshold = value;
        }

        /// <summary>
        /// Angular speed (rad/s) below which a body counts as resting.
        /// </summary>
        public double SleepAngularVelocityThreshold
        {
            get => _sleep.SleepAngularVelocityThreshold;
            set => _sleep.SleepAngularVelocityThreshold = value;
        }

        /// <summary>
        /// Time (s) a whole island has to rest before it is put to sleep.
        /// </summary>
        public double SleepTimeThreshold
        {
            get => _sleep.SleepTimeThreshold;
            set => _sleep.SleepTimeThreshold = value;
        }

        /// <summary>
        /// Edge length of the broad-phase grid cells. Should be in the order of the
        /// size of a typical body.
        /// </summary>
        public double BroadPhaseCellSize
        {
            get => _broadPhaseCellSize;
            set
            {
                if (!(value > 0) || double.IsInfinity(value))
                    throw new ArgumentOutOfRangeException(nameof(value), "Cell size must be positive and finite.");
                _broadPhaseCellSize = value;
            }
        }

        /// <summary>
        /// Gets or sets whether physics simulation is paused.
        /// </summary>
        public bool IsPaused { get; set; }

        /// <summary>
        /// Gets or sets the velocity threshold for collision restitution.
        /// </summary>
        public double RestitutionThreshold { get; set; } = 1.0;

        /// <summary>
        /// Gets or sets the world bounds. Bodies outside may be deactivated.
        /// </summary>
        public AABB? WorldBounds { get; set; }

        /// <summary>
        /// Gets or sets whether to use continuous collision detection. When enabled,
        /// fast sphere and box bodies are stopped at their first time of impact instead
        /// of tunnelling through thin obstacles.
        /// </summary>
        public bool UseCCD { get; set; } = false;

        /// <summary>
        /// Gets the total simulation time.
        /// </summary>
        public double TotalTime { get; private set; }

        /// <summary>
        /// Gets the number of physics steps taken.
        /// </summary>
        public long StepCount { get; private set; }

        /// <summary>
        /// Event raised when a collision occurs.
        /// </summary>
        public event Action<CollisionInfo>? OnCollision;

        /// <summary>
        /// Event raised before each physics step.
        /// </summary>
        public event Action<double>? OnPreStep;

        /// <summary>
        /// Event raised after each physics step.
        /// </summary>
        public event Action<double>? OnPostStep;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new physics world with default settings.
        /// </summary>
        public PhysicsWorld()
        {
            _bodies = new List<IPhysicsBody>();
            _globalForces = new List<IForce>();
            _particleSystems = new List<ParticleSystem>();
            _collisions = new List<CollisionInfo>();
            _spatialGrid = new Dictionary<long, List<int>>();
        }

        /// <summary>
        /// Creates a new physics world with custom gravity.
        /// </summary>
        public PhysicsWorld(Vector3D gravity) : this()
        {
            Gravity = gravity;
        }

        #endregion

        #region Body Management

        /// <summary>
        /// Adds a body to the world.
        /// </summary>
        public void AddBody(IPhysicsBody body)
        {
            _bodiesLock.EnterWriteLock();
            try
            {
                if (!_bodies.Contains(body))
                {
                    _bodies.Add(body);
                    _bodySnapshot = null;
                }
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        /// <summary>
        /// Removes a body from the world.
        /// </summary>
        public bool RemoveBody(IPhysicsBody body)
        {
            bool removed;
            _bodiesLock.EnterWriteLock();
            try
            {
                removed = _bodies.Remove(body);
                if (removed)
                    _bodySnapshot = null;
            }
            finally { _bodiesLock.ExitWriteLock(); }

            // Whatever was resting on the removed body has lost its support.
            if (removed)
                WakeAllBodies();

            return removed;
        }

        /// <summary>
        /// Gets a body by its ID.
        /// </summary>
        public IPhysicsBody? GetBody(string id)
        {
            _bodiesLock.EnterReadLock();
            try { return _bodies.Find(b => b.Id == id); }
            finally { _bodiesLock.ExitReadLock(); }
        }

        /// <summary>
        /// Clears all bodies from the world.
        /// </summary>
        public void ClearBodies()
        {
            _bodiesLock.EnterWriteLock();
            try
            {
                _bodies.Clear();
                _bodySnapshot = null;
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        #endregion

        #region Force Management

        public void AddForce(IForce force)
        {
            if (!_globalForces.Contains(force))
                _globalForces.Add(force);
        }

        public bool RemoveForce(IForce force) => _globalForces.Remove(force);
        public void ClearForces() => _globalForces.Clear();

        #endregion

        #region Particle System Management

        public void AddParticleSystem(ParticleSystem system)
        {
            if (!_particleSystems.Contains(system))
                _particleSystems.Add(system);
        }

        public bool RemoveParticleSystem(ParticleSystem system) => _particleSystems.Remove(system);

        #endregion

        #region Simulation (Multi-threaded)

        /// <summary>
        /// Updates the physics simulation with multi-threading.
        /// </summary>
        public void Update(double deltaTime)
        {
            if (IsPaused)
                return;

            _accumulator += deltaTime;

            double maxAccumulator = FixedTimeStep * MaxSubSteps;
            if (_accumulator > maxAccumulator)
                _accumulator = maxAccumulator;

            while (_accumulator >= FixedTimeStep)
            {
                Step(FixedTimeStep);
                _accumulator -= FixedTimeStep;
            }

            // Update particle systems in parallel
            Parallel.ForEach(_particleSystems, system => system.Update(deltaTime));
        }

        /// <summary>
        /// Performs a single physics step.
        /// <para>
        /// Pipeline: forces -> velocities -> collision detection -> contact solver ->
        /// positions (with optional continuous collision detection) -> penetration
        /// recovery -> sleeping.
        /// </para>
        /// </summary>
        public void Step(double dt)
        {
            OnPreStep?.Invoke(dt);

            _collisions.Clear();

            var bodies = GetBodySnapshot();

            // 1. Forces and velocities. Positions are integrated only after the
            //    contact solver has removed approaching velocities.
            IntegrateVelocities(bodies, dt);

            // 2. Collision detection at the current positions.
            BroadPhase(bodies);
            NarrowPhase(bodies);

            // 3. Wake sleeping bodies that are being disturbed, then solve contacts.
            _sleep.BuildIslands(bodies, _collisions);
            _sleep.WakeDisturbedIslands(bodies, _collisions);

            _solver.VelocityIterations = SolverIterations;
            _solver.PositionIterations = PositionIterations;
            _solver.RestitutionThreshold = RestitutionThreshold;
            _solver.SolveVelocities(_collisions);

            // 4. Positions.
            IntegratePositions(bodies, dt);
            _solver.SolvePositions();

            // 5. Sleeping.
            if (AllowSleeping)
                _sleep.UpdateSleeping(bodies, dt);

            if (OnCollision != null)
            {
                foreach (var collision in _collisions)
                    OnCollision.Invoke(collision);
            }

            TotalTime += dt;
            StepCount++;

            OnPostStep?.Invoke(dt);
        }

        private IPhysicsBody[] GetBodySnapshot()
        {
            var snapshot = _bodySnapshot;
            if (snapshot != null)
                return snapshot;

            _bodiesLock.EnterReadLock();
            try
            {
                snapshot = _bodies.ToArray();
                _bodySnapshot = snapshot;
                return snapshot;
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        private void WakeAllBodies()
        {
            foreach (var body in GetBodySnapshot())
            {
                if (body.IsSleeping)
                    body.WakeUp();
            }
        }

        private void IntegrateVelocities(IPhysicsBody[] bodies, double dt)
        {
            var gravity = Gravity;
            var forces = _globalForces.Count > 0 ? _globalForces.ToArray() : Array.Empty<IForce>();

            void Process(IPhysicsBody body)
            {
                if (!body.IsActive)
                    return;

                if (body.BodyType == BodyType.Dynamic && !body.IsSleeping)
                {
                    body.ApplyForce(gravity * body.Mass);

                    foreach (var force in forces)
                    {
                        if (force.Enabled)
                            body.ApplyForce(force.Calculate(body.Position, body.Velocity, body.Mass));
                    }
                }

                if (body is RigidBody rb)
                    rb.IntegrateVelocity(dt);
                else
                    body.Integrate(dt);
            }

            if (bodies.Length >= ParallelThreshold)
                Parallel.ForEach(bodies, Process);
            else
                foreach (var body in bodies) Process(body);
        }

        private void IntegratePositions(IPhysicsBody[] bodies, double dt)
        {
            bool useCcd = UseCCD;

            // Sequential: the continuous collision query reads other bodies' positions.
            foreach (var body in bodies)
            {
                if (body is not RigidBody rb || !rb.IsActive || rb.IsSleeping || rb.BodyType == BodyType.Static)
                    continue;

                double stepTime = dt;
                if (useCcd && rb.BodyType == BodyType.Dynamic && ContinuousCollision.NeedsSweep(rb, dt, out var swept))
                    stepTime = ContinuousCollision.TimeOfImpact(rb, CollectSweepCandidates(bodies, swept), dt, swept);

                rb.IntegratePosition(stepTime);
            }
        }

        #endregion

        #region Broad Phase

        /// <summary>
        /// Bodies near a swept volume, taken from the broad-phase grid of this step.
        /// The result may contain duplicates.
        /// </summary>
        private IReadOnlyList<IPhysicsBody> CollectSweepCandidates(IPhysicsBody[] bodies, AABB swept)
        {
            // Bodies integrated earlier in this step have already moved a little.
            swept = swept.Expanded(0.5 * BroadPhaseCellSize);

            double inv = 1.0 / BroadPhaseCellSize;
            double minX = Math.Floor(swept.Min.X * inv), maxX = Math.Floor(swept.Max.X * inv);
            double minY = Math.Floor(swept.Min.Y * inv), maxY = Math.Floor(swept.Max.Y * inv);
            double minZ = Math.Floor(swept.Min.Z * inv), maxZ = Math.Floor(swept.Max.Z * inv);
            double cells = (maxX - minX + 1) * (maxY - minY + 1) * (maxZ - minZ + 1);

            // A sweep across a huge volume is cheaper to test against everything.
            if (!(cells <= 4096) || Math.Abs(minX) > CellCoordLimit || Math.Abs(maxX) > CellCoordLimit ||
                Math.Abs(minY) > CellCoordLimit || Math.Abs(maxY) > CellCoordLimit ||
                Math.Abs(minZ) > CellCoordLimit || Math.Abs(maxZ) > CellCoordLimit)
                return bodies;

            _sweepCandidates.Clear();
            foreach (int large in _largeBodies)
                _sweepCandidates.Add(bodies[large]);

            for (int cx = (int)minX; cx <= (int)maxX; cx++)
            {
                for (int cy = (int)minY; cy <= (int)maxY; cy++)
                {
                    for (int cz = (int)minZ; cz <= (int)maxZ; cz++)
                    {
                        if (!_spatialGrid.TryGetValue(GetCellKeyFromCoords(cx, cy, cz), out var list))
                            continue;

                        foreach (int index in list)
                            _sweepCandidates.Add(bodies[index]);
                    }
                }
            }

            return _sweepCandidates;
        }

        private void BroadPhase(IPhysicsBody[] bodies)
        {
            foreach (var list in _spatialGrid.Values)
            {
                list.Clear();
                _cellListPool.Push(list);
            }
            _spatialGrid.Clear();
            _largeBodies.Clear();
            _pairs.Clear();

            int count = bodies.Length;
            if (_cellMin.Length < count)
            {
                _cellMin = new (int, int, int)[Math.Max(count, _cellMin.Length * 2)];
                _cellMax = new (int, int, int)[_cellMin.Length];
            }

            double inv = 1.0 / BroadPhaseCellSize;

            for (int i = 0; i < count; i++)
            {
                var body = bodies[i];
                if (!body.IsActive)
                    continue;

                var box = body.BoundingBox;
                double minX = Math.Floor(box.Min.X * inv), maxX = Math.Floor(box.Max.X * inv);
                double minY = Math.Floor(box.Min.Y * inv), maxY = Math.Floor(box.Max.Y * inv);
                double minZ = Math.Floor(box.Min.Z * inv), maxZ = Math.Floor(box.Max.Z * inv);

                double cells = (maxX - minX + 1) * (maxY - minY + 1) * (maxZ - minZ + 1);

                // Bodies that span many cells (ground planes, walls) or whose bounds are
                // not finite are tested against everything instead of flooding the grid.
                if (!(cells <= MaxCellsPerBody) || Math.Abs(minX) > CellCoordLimit || Math.Abs(maxX) > CellCoordLimit ||
                    Math.Abs(minY) > CellCoordLimit || Math.Abs(maxY) > CellCoordLimit ||
                    Math.Abs(minZ) > CellCoordLimit || Math.Abs(maxZ) > CellCoordLimit)
                {
                    if (double.IsNaN(cells))
                        continue;
                    _largeBodies.Add(i);
                    continue;
                }

                var min = ((int)minX, (int)minY, (int)minZ);
                var max = ((int)maxX, (int)maxY, (int)maxZ);
                _cellMin[i] = min;
                _cellMax[i] = max;

                for (int cx = min.Item1; cx <= max.Item1; cx++)
                {
                    for (int cy = min.Item2; cy <= max.Item2; cy++)
                    {
                        for (int cz = min.Item3; cz <= max.Item3; cz++)
                        {
                            long key = GetCellKeyFromCoords(cx, cy, cz);
                            if (!_spatialGrid.TryGetValue(key, out var list))
                            {
                                list = _cellListPool.Count > 0 ? _cellListPool.Pop() : new List<int>(8);
                                _spatialGrid[key] = list;
                            }
                            list.Add(i);
                        }
                    }
                }
            }

            // Pairs sharing a cell. Two bodies can share several cells, so a pair is
            // reported only from the first cell their cell ranges have in common.
            foreach (var cell in _spatialGrid)
            {
                var indices = cell.Value;
                if (indices.Count < 2)
                    continue;

                long cellKey = cell.Key;

                for (int a = 0; a < indices.Count; a++)
                {
                    int i = indices[a];
                    for (int b = a + 1; b < indices.Count; b++)
                    {
                        int j = indices[b];

                        var minI = _cellMin[i];
                        var minJ = _cellMin[j];
                        long firstShared = GetCellKeyFromCoords(
                            Math.Max(minI.Item1, minJ.Item1),
                            Math.Max(minI.Item2, minJ.Item2),
                            Math.Max(minI.Item3, minJ.Item3));
                        if (firstShared != cellKey)
                            continue;

                        AddPairIfOverlapping(bodies, i, j);
                    }
                }
            }

            // Large bodies against everything else.
            for (int a = 0; a < _largeBodies.Count; a++)
            {
                int large = _largeBodies[a];
                for (int k = 0; k < count; k++)
                {
                    if (k == large)
                        continue;

                    // Large-large pairs are visited from both sides: keep one.
                    if (k < large && _largeBodies.Contains(k))
                        continue;

                    AddPairIfOverlapping(bodies, large, k);
                }
            }
        }

        private void AddPairIfOverlapping(IPhysicsBody[] bodies, int i, int j)
        {
            var bodyA = bodies[i];
            var bodyB = bodies[j];

            if (!bodyA.IsActive || !bodyB.IsActive)
                return;

            // At least one side must be able to react.
            bool movingA = bodyA.BodyType == BodyType.Dynamic && !bodyA.IsSleeping;
            bool movingB = bodyB.BodyType == BodyType.Dynamic && !bodyB.IsSleeping;
            if (!movingA && !movingB)
            {
                // A moving kinematic body must still be able to wake a sleeper.
                bool sleeperAndMover =
                    (bodyA.BodyType == BodyType.Dynamic && SleepManager.IsMovingKinematic(bodyB)) ||
                    (bodyB.BodyType == BodyType.Dynamic && SleepManager.IsMovingKinematic(bodyA));
                if (!sleeperAndMover)
                    return;
            }

            if (!bodyA.BoundingBox.Intersects(bodyB.BoundingBox))
                return;

            _pairs.Add(i < j ? (i, j) : (j, i));
        }

        private static long GetCellKeyFromCoords(int cx, int cy, int cz)
        {
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (long)(cz & 0x1FFFFF);
        }

        private void NarrowPhase(IPhysicsBody[] bodies)
        {
            int pairCount = _pairs.Count;
            if (pairCount == 0)
                return;

            if (pairCount >= ParallelThreshold)
            {
                // Results are written by pair index so the collision order (and with it
                // the simulation) does not depend on thread scheduling.
                if (_narrowPhaseResults.Length < pairCount)
                    _narrowPhaseResults = new CollisionInfo[Math.Max(pairCount, _narrowPhaseResults.Length * 2)];

                var results = _narrowPhaseResults;
                var pairs = _pairs;
                Parallel.For(0, pairCount, p =>
                {
                    var (i, j) = pairs[p];
                    results[p] = CollisionDetector.Detect(bodies[i], bodies[j]);
                });

                for (int p = 0; p < pairCount; p++)
                {
                    if (results[p].HasCollision)
                        _collisions.Add(results[p]);
                    results[p] = default;
                }
            }
            else
            {
                for (int p = 0; p < pairCount; p++)
                {
                    var (i, j) = _pairs[p];
                    var collision = CollisionDetector.Detect(bodies[i], bodies[j]);
                    if (collision.HasCollision)
                        _collisions.Add(collision);
                }
            }
        }

        #endregion

        #region Queries

        /// <summary>
        /// Casts a ray through the world.
        /// </summary>
        public RaycastHit? Raycast(Vector3D origin, Vector3D direction, double maxDistance = 1000)
        {
            RaycastHit? closest = null;
            double closestDist = maxDistance;

            _bodiesLock.EnterReadLock();
            try
            {
                foreach (var body in _bodies)
                {
                    if (!body.IsActive) continue;

                    if (body is RigidBody rb)
                    {
                        if (CollisionDetector.Raycast(origin, direction, rb, closestDist,
                            out double dist, out Vector3D point, out Vector3D normal))
                        {
                            if (dist < closestDist)
                            {
                                closestDist = dist;
                                closest = new RaycastHit
                                {
                                    Body = body,
                                    Distance = dist,
                                    Point = point,
                                    Normal = normal
                                };
                            }
                        }
                    }
                }
            }
            finally { _bodiesLock.ExitReadLock(); }

            return closest;
        }

        /// <summary>
        /// Gets all bodies within a sphere (parallelized).
        /// </summary>
        public List<IPhysicsBody> QuerySphere(Vector3D center, double radius)
        {
            var results = new ConcurrentBag<IPhysicsBody>();
            double radiusSq = radius * radius;

            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (!body.IsActive) return;

                    if (Vector3D.DistanceSquared(center, body.Position) <= radiusSq)
                    {
                        results.Add(body);
                    }
                });
            }
            finally { _bodiesLock.ExitReadLock(); }

            return new List<IPhysicsBody>(results);
        }

        /// <summary>
        /// Gets all bodies within an AABB (parallelized).
        /// </summary>
        public List<IPhysicsBody> QueryAABB(AABB bounds)
        {
            var results = new ConcurrentBag<IPhysicsBody>();

            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (!body.IsActive) return;

                    if (bounds.Intersects(body.BoundingBox))
                    {
                        results.Add(body);
                    }
                });
            }
            finally { _bodiesLock.ExitReadLock(); }

            return new List<IPhysicsBody>(results);
        }

        #endregion

        #region Utility

        /// <summary>
        /// Resets the simulation.
        /// </summary>
        public void Reset()
        {
            TotalTime = 0;
            StepCount = 0;
            _accumulator = 0;
            _collisions.Clear();
            _solver.ClearCache();

            _bodiesLock.EnterReadLock();
            try
            {
                foreach (var body in _bodies)
                {
                    body.WakeUp();
                }
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        /// <summary>
        /// Creates a static ground plane.
        /// </summary>
        public RigidBody CreateGround(double width = 100, double depth = 100, double height = 1)
        {
            var ground = RigidBody.CreateStaticBox(
                new Vector3D(0, -height / 2, 0),
                new Vector3D(width / 2, height / 2, depth / 2)
            );
            AddBody(ground);
            return ground;
        }

        /// <summary>
        /// Creates walls around an area.
        /// </summary>
        public void CreateWalls(double width, double height, double depth, double thickness = 1)
        {
            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(-width / 2 - thickness / 2, height / 2, 0),
                new Vector3D(thickness / 2, height / 2, depth / 2)));

            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(width / 2 + thickness / 2, height / 2, 0),
                new Vector3D(thickness / 2, height / 2, depth / 2)));

            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(0, height / 2, depth / 2 + thickness / 2),
                new Vector3D(width / 2, height / 2, thickness / 2)));

            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(0, height / 2, -depth / 2 - thickness / 2),
                new Vector3D(width / 2, height / 2, thickness / 2)));
        }

        #endregion
    }

    /// <summary>
    /// Information about a raycast hit.
    /// </summary>
    public struct RaycastHit
    {
        public IPhysicsBody Body;
        public double Distance;
        public Vector3D Point;
        public Vector3D Normal;
    }
}
