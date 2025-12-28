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
        private readonly ConcurrentBag<CollisionInfo> _collisionBag;
        private readonly List<CollisionInfo> _collisions;
        private double _accumulator;

        // Spatial hashing for broad phase
        private readonly Dictionary<long, List<int>> _spatialGrid;
        private const float CellSize = 2.0f;
        private const float InvCellSize = 1f / CellSize;

        // Thread synchronization
        private readonly ReaderWriterLockSlim _bodiesLock = new();

        #endregion

        #region Properties

        /// <summary>
        /// Gets all bodies in the world.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> Bodies
        {
            get
            {
                _bodiesLock.EnterReadLock();
                try { return _bodies.ToArray(); }
                finally { _bodiesLock.ExitReadLock(); }
            }
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
        public Vector3D Gravity { get; set; } = new(0, -PhysicsConstants.EarthGravity, 0);

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
        /// Gets or sets whether to use continuous collision detection.
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
            _collisionBag = new ConcurrentBag<CollisionInfo>();
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
                    _bodies.Add(body);
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        /// <summary>
        /// Removes a body from the world.
        /// </summary>
        public bool RemoveBody(IPhysicsBody body)
        {
            _bodiesLock.EnterWriteLock();
            try { return _bodies.Remove(body); }
            finally { _bodiesLock.ExitWriteLock(); }
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
            try { _bodies.Clear(); }
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
        /// Performs a single physics step with parallel processing.
        /// </summary>
        public void Step(double dt)
        {
            OnPreStep?.Invoke(dt);

            _collisions.Clear();

            // Apply forces in parallel
            ApplyForcesParallel(dt);

            // Integrate in parallel
            IntegrateBodiesParallel(dt);

            // Build spatial grid and detect collisions
            BuildSpatialGrid();
            DetectCollisionsParallel();

            // Resolve collisions (sequential for determinism)
            ResolveCollisions();

            TotalTime += dt;
            StepCount++;

            OnPostStep?.Invoke(dt);
        }

        private void ApplyForcesParallel(double dt)
        {
            _bodiesLock.EnterReadLock();
            try
            {
                var gravity = Gravity;
                var forces = _globalForces.ToArray();

                Parallel.ForEach(_bodies, body =>
                {
                    if (body.BodyType != BodyType.Dynamic || !body.IsActive || body.IsSleeping)
                        return;

                    body.ApplyForce(gravity * body.Mass);

                    foreach (var force in forces)
                    {
                        if (force.Enabled)
                        {
                            var f = force.Calculate(body.Position, body.Velocity, body.Mass);
                            body.ApplyForce(f);
                        }
                    }
                });
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        private void IntegrateBodiesParallel(double dt)
        {
            _bodiesLock.EnterReadLock();
            try
            {
                var bounds = WorldBounds;

                Parallel.ForEach(_bodies, body =>
                {
                    if (!body.IsActive) return;

                    body.Integrate(dt);

                    // Check world bounds (optional deactivation)
                    // if (bounds.HasValue && !bounds.Value.Contains(body.Position))
                    //     body.IsActive = false;
                });
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        private void BuildSpatialGrid()
        {
            _spatialGrid.Clear();

            _bodiesLock.EnterReadLock();
            try
            {
                for (int i = 0; i < _bodies.Count; i++)
                {
                    var body = _bodies[i];
                    if (!body.IsActive) continue;

                    var pos = body.Position;
                    long cellKey = GetCellKey((float)pos.X, (float)pos.Y, (float)pos.Z);

                    if (!_spatialGrid.TryGetValue(cellKey, out var list))
                    {
                        list = new List<int>(8);
                        _spatialGrid[cellKey] = list;
                    }
                    list.Add(i);
                }
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        private static long GetCellKey(float x, float y, float z)
        {
            int cx = (int)MathF.Floor(x * InvCellSize);
            int cy = (int)MathF.Floor(y * InvCellSize);
            int cz = (int)MathF.Floor(z * InvCellSize);
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (cz & 0x1FFFFF);
        }

        private static long GetCellKeyFromCoords(int cx, int cy, int cz)
        {
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (cz & 0x1FFFFF);
        }

        private void DetectCollisionsParallel()
        {
            // Clear the concurrent bag
            while (_collisionBag.TryTake(out _)) { }

            _bodiesLock.EnterReadLock();
            try
            {
                var cells = _spatialGrid.ToArray();
                var bodiesArray = _bodies.ToArray();

                Parallel.ForEach(cells, cell =>
                {
                    var indices = cell.Value;
                    long cellKey = cell.Key;

                    // Unpack cell coordinates
                    int cx = (int)((cellKey >> 42) & 0x1FFFFF);
                    int cy = (int)((cellKey >> 21) & 0x1FFFFF);
                    int cz = (int)(cellKey & 0x1FFFFF);
                    if (cx >= 0x100000) cx -= 0x200000;
                    if (cy >= 0x100000) cy -= 0x200000;
                    if (cz >= 0x100000) cz -= 0x200000;

                    // Check within cell
                    for (int i = 0; i < indices.Count; i++)
                    {
                        for (int j = i + 1; j < indices.Count; j++)
                        {
                            CheckCollisionPair(bodiesArray, indices[i], indices[j]);
                        }
                    }

                    // Check neighboring cells (positive direction only to avoid duplicates)
                    for (int dx = 0; dx <= 1; dx++)
                    {
                        for (int dy = 0; dy <= 1; dy++)
                        {
                            for (int dz = 0; dz <= 1; dz++)
                            {
                                if (dx == 0 && dy == 0 && dz == 0) continue;

                                long neighborKey = GetCellKeyFromCoords(cx + dx, cy + dy, cz + dz);
                                if (!_spatialGrid.TryGetValue(neighborKey, out var neighborIndices)) continue;

                                foreach (var i in indices)
                                {
                                    foreach (var j in neighborIndices)
                                    {
                                        CheckCollisionPair(bodiesArray, i, j);
                                    }
                                }
                            }
                        }
                    }
                });

                // Transfer from concurrent bag to list
                while (_collisionBag.TryTake(out var collision))
                {
                    _collisions.Add(collision);
                }
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        private void CheckCollisionPair(IPhysicsBody[] bodies, int i, int j)
        {
            var bodyA = bodies[i];
            var bodyB = bodies[j];

            if (!bodyA.IsActive || !bodyB.IsActive) return;
            if (bodyA.BodyType == BodyType.Static && bodyB.BodyType == BodyType.Static) return;

            var collision = CollisionDetector.Detect(bodyA, bodyB);
            if (collision.HasCollision)
            {
                _collisionBag.Add(collision);
            }
        }

        private void ResolveCollisions()
        {
            foreach (var collision in _collisions)
            {
                CollisionResolver.Resolve(collision, RestitutionThreshold);
                OnCollision?.Invoke(collision);
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
