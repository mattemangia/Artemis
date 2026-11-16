using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Runtime.CompilerServices;
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
    /// High-performance physics world optimized for real-time applications.
    /// Features:
    /// - Spatial hashing for O(1) average-case broad-phase collision detection
    /// - Island solver for parallel constraint solving
    /// - Object pooling to reduce GC pressure
    /// - Warm starting for faster convergence
    /// - Multi-threaded processing
    /// - Continuous collision detection (CCD)
    /// </summary>
    public class OptimizedPhysicsWorld
    {
        #region Constants

        private const double SleepVelocityThreshold = 0.1;
        private const double SleepAngularThreshold = 0.25;
        private const double SleepTimeRequired = 0.5;
        private const int DefaultSolverIterations = 8;
        private const int DefaultPositionIterations = 3;

        #endregion

        #region Fields

        private readonly List<IPhysicsBody> _bodies;
        private readonly List<IPhysicsBody> _dynamicBodies;
        private readonly List<IForce> _globalForces;
        private readonly List<ParticleSystem> _particleSystems;
        private readonly List<CollisionInfo> _collisions;
        private readonly List<(IPhysicsBody A, IPhysicsBody B)> _broadPhasePairs;

        // Optimizations
        private readonly SpatialHash _spatialHash;
        private readonly IslandSolver _islandSolver;
        private readonly ListPool<CollisionInfo> _collisionListPool;

        // Contact solving (with warm starting) and island sleeping
        private readonly ContactSolver _solver = new();
        private readonly SleepManager _sleep = new()
        {
            SleepVelocityThreshold = SleepVelocityThreshold,
            SleepAngularVelocityThreshold = SleepAngularThreshold,
            SleepTimeThreshold = SleepTimeRequired
        };

        // Threading
        private readonly int _threadCount;
        private readonly ParallelOptions _parallelOptions;

        // Timing
        private double _accumulator;
        private readonly Stopwatch _profileTimer;

        #endregion

        #region Properties

        /// <summary>Gets all bodies in the world.</summary>
        public IReadOnlyList<IPhysicsBody> Bodies => _bodies;

        /// <summary>Gets all dynamic bodies (optimized access).</summary>
        public IReadOnlyList<IPhysicsBody> DynamicBodies => _dynamicBodies;

        /// <summary>Gets all global forces.</summary>
        public IReadOnlyList<IForce> GlobalForces => _globalForces;

        /// <summary>Gets all particle systems.</summary>
        public IReadOnlyList<ParticleSystem> ParticleSystems => _particleSystems;

        /// <summary>Gets the collisions from the last update.</summary>
        public IReadOnlyList<CollisionInfo> LastCollisions => _collisions;

        /// <summary>Gets or sets the fixed time step.</summary>
        public double FixedTimeStep { get; set; } = 1.0 / 60.0;

        /// <summary>Gets or sets the maximum number of sub-steps per frame.</summary>
        public int MaxSubSteps { get; set; } = 8;

        /// <summary>Gets or sets the gravity vector.</summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>Gets or sets whether simulation is paused.</summary>
        public bool IsPaused { get; set; }

        /// <summary>Gets or sets the number of solver iterations.</summary>
        public int SolverIterations { get; set; } = DefaultSolverIterations;

        /// <summary>Gets or sets the number of position correction iterations.</summary>
        public int PositionIterations { get; set; } = DefaultPositionIterations;

        /// <summary>Gets or sets the velocity threshold for collision restitution.</summary>
        public double RestitutionThreshold { get; set; } = 1.0;

        /// <summary>Gets or sets whether to use continuous collision detection.</summary>
        public bool UseCCD { get; set; } = true;

        /// <summary>Gets or sets whether to use warm starting.</summary>
        public bool UseWarmStarting { get; set; } = true;

        /// <summary>Gets or sets whether to use multi-threading.</summary>
        public bool UseMultiThreading { get; set; } = true;

        /// <summary>Gets or sets the world bounds.</summary>
        public AABB? WorldBounds { get; set; }

        /// <summary>Gets the total simulation time.</summary>
        public double TotalTime { get; private set; }

        /// <summary>Gets the number of physics steps taken.</summary>
        public long StepCount { get; private set; }

        /// <summary>Gets the island solver for diagnostics.</summary>
        public IslandSolver IslandSolver => _islandSolver;

        /// <summary>Gets the spatial hash for diagnostics.</summary>
        public SpatialHash SpatialHash => _spatialHash;

        // Performance metrics
        /// <summary>Broad phase time in milliseconds.</summary>
        public double BroadPhaseTimeMs { get; private set; }

        /// <summary>Narrow phase time in milliseconds.</summary>
        public double NarrowPhaseTimeMs { get; private set; }

        /// <summary>Solver time in milliseconds.</summary>
        public double SolverTimeMs { get; private set; }

        /// <summary>Integration time in milliseconds.</summary>
        public double IntegrationTimeMs { get; private set; }

        /// <summary>Event raised when a collision occurs.</summary>
        public event Action<CollisionInfo>? OnCollision;

        /// <summary>Event raised before each physics step.</summary>
        public event Action<double>? OnPreStep;

        /// <summary>Event raised after each physics step.</summary>
        public event Action<double>? OnPostStep;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new optimized physics world.
        /// </summary>
        /// <param name="spatialCellSize">Cell size for spatial hash (should be >= largest object).</param>
        /// <param name="threadCount">Number of threads (0 = auto-detect).</param>
        public OptimizedPhysicsWorld(double spatialCellSize = 2.0, int threadCount = 0)
        {
            _bodies = new List<IPhysicsBody>(1024);
            _dynamicBodies = new List<IPhysicsBody>(512);
            _globalForces = new List<IForce>(8);
            _particleSystems = new List<ParticleSystem>(16);
            _collisions = new List<CollisionInfo>(256);
            _broadPhasePairs = new List<(IPhysicsBody, IPhysicsBody)>(512);

            _spatialHash = new SpatialHash(spatialCellSize);
            _islandSolver = new IslandSolver(10000);
            _collisionListPool = new ListPool<CollisionInfo>(64);

            _threadCount = threadCount > 0 ? threadCount : Environment.ProcessorCount;
            _parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = _threadCount
            };

            _profileTimer = new Stopwatch();
        }

        /// <summary>
        /// Creates a new optimized physics world with custom gravity.
        /// </summary>
        public OptimizedPhysicsWorld(Vector3D gravity, double spatialCellSize = 2.0)
            : this(spatialCellSize)
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
            if (!_bodies.Contains(body))
            {
                _bodies.Add(body);
                if (body.BodyType == BodyType.Dynamic)
                {
                    _dynamicBodies.Add(body);
                }
            }
        }

        /// <summary>
        /// Removes a body from the world.
        /// </summary>
        public bool RemoveBody(IPhysicsBody body)
        {
            _dynamicBodies.Remove(body);
            bool removed = _bodies.Remove(body);

            // Whatever was resting on the removed body has lost its support.
            if (removed)
            {
                foreach (var other in _bodies)
                {
                    if (other.IsSleeping)
                        other.WakeUp();
                }
            }

            return removed;
        }

        /// <summary>
        /// Gets a body by its ID.
        /// </summary>
        public IPhysicsBody? GetBody(string id)
        {
            return _bodies.Find(b => b.Id == id);
        }

        /// <summary>
        /// Clears all bodies from the world.
        /// </summary>
        public void ClearBodies()
        {
            _bodies.Clear();
            _dynamicBodies.Clear();
            _spatialHash.Clear();
            _solver.ClearCache();
        }

        /// <summary>
        /// Adds a global force.
        /// </summary>
        public void AddForce(IForce force)
        {
            if (!_globalForces.Contains(force))
            {
                _globalForces.Add(force);
            }
        }

        /// <summary>
        /// Removes a global force.
        /// </summary>
        public bool RemoveForce(IForce force)
        {
            return _globalForces.Remove(force);
        }

        /// <summary>
        /// Adds a particle system.
        /// </summary>
        public void AddParticleSystem(ParticleSystem system)
        {
            if (!_particleSystems.Contains(system))
            {
                _particleSystems.Add(system);
            }
        }

        /// <summary>
        /// Removes a particle system.
        /// </summary>
        public bool RemoveParticleSystem(ParticleSystem system)
        {
            return _particleSystems.Remove(system);
        }

        #endregion

        #region Simulation

        /// <summary>
        /// Updates the physics simulation.
        /// </summary>
        /// <param name="deltaTime">Time since last update in seconds.</param>
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

            // Update particle systems
            foreach (var system in _particleSystems)
            {
                system.Update(deltaTime);
            }
        }

        /// <summary>
        /// Performs a single physics step.
        /// </summary>
        public void Step(double dt)
        {
            OnPreStep?.Invoke(dt);

            _collisions.Clear();

            // 1. Apply forces
            _profileTimer.Restart();
            ApplyForces(dt);

            // 2. Integrate velocities
            IntegrateVelocities(dt);

            IntegrationTimeMs = _profileTimer.Elapsed.TotalMilliseconds;

            // 3. Broad phase collision detection
            _profileTimer.Restart();
            BroadPhase();
            BroadPhaseTimeMs = _profileTimer.Elapsed.TotalMilliseconds;

            // 4. Narrow phase collision detection
            _profileTimer.Restart();
            NarrowPhase();
            NarrowPhaseTimeMs = _profileTimer.Elapsed.TotalMilliseconds;

            // 5. Build islands and wake the ones that are being disturbed
            _islandSolver.BuildIslands(_bodies, _collisions);
            _sleep.BuildIslands(_bodies, _collisions);
            _sleep.WakeDisturbedIslands(_bodies, _collisions);

            // 6. Solve constraints
            _profileTimer.Restart();
            SolveConstraints(dt);
            SolverTimeMs = _profileTimer.Elapsed.TotalMilliseconds;

            // 7. Integrate positions
            IntegratePositions(dt);

            // 8. Update sleeping
            _sleep.UpdateSleeping(_bodies, dt);

            // 9. Fire collision events
            foreach (var collision in _collisions)
            {
                OnCollision?.Invoke(collision);
            }

            TotalTime += dt;
            StepCount++;

            OnPostStep?.Invoke(dt);
        }

        #endregion

        #region Physics Pipeline

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void ApplyForces(double dt)
        {
            if (UseMultiThreading && _bodies.Count > 100)
            {
                Parallel.ForEach(_bodies, _parallelOptions, body =>
                {
                    ApplyForcesToBody(body);
                });
            }
            else
            {
                foreach (var body in _bodies)
                {
                    ApplyForcesToBody(body);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void ApplyForcesToBody(IPhysicsBody body)
        {
            // The body type can change after a body was added, so it is checked here
            // rather than relying on the list the body was sorted into.
            if (body.BodyType != BodyType.Dynamic || !body.IsActive || body.IsSleeping)
                return;

            // Apply gravity
            body.ApplyForce(Gravity * body.Mass);

            // Apply global forces
            foreach (var force in _globalForces)
            {
                if (force.Enabled)
                {
                    var f = force.Calculate(body.Position, body.Velocity, body.Mass);
                    body.ApplyForce(f);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void IntegrateVelocities(double dt)
        {
            if (UseMultiThreading && _bodies.Count > 100)
            {
                Parallel.ForEach(_bodies, _parallelOptions, body =>
                {
                    if (body.IsActive && !body.IsSleeping)
                    {
                        IntegrateVelocity(body, dt);
                    }
                });
            }
            else
            {
                foreach (var body in _bodies)
                {
                    if (body.IsActive && !body.IsSleeping)
                    {
                        IntegrateVelocity(body, dt);
                    }
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void IntegrateVelocity(IPhysicsBody body, double dt)
        {
            // Forces -> velocities only. Positions follow after the contact solver.
            if (body is RigidBody rb)
                rb.IntegrateVelocity(dt);
            else
                body.Integrate(dt);
        }

        private void BroadPhase()
        {
            // Rebuild spatial hash
            _spatialHash.Rebuild(_bodies);

            // Get potential collision pairs
            _spatialHash.GetPotentialPairs(_broadPhasePairs);
        }

        private void NarrowPhase()
        {
            if (UseMultiThreading && _broadPhasePairs.Count > 50)
            {
                var localCollisions = new ConcurrentBag<CollisionInfo>();

                Parallel.ForEach(_broadPhasePairs, _parallelOptions, pair =>
                {
                    var collision = CollisionDetector.Detect(pair.A, pair.B);
                    if (collision.HasCollision)
                    {
                        localCollisions.Add(collision);
                    }
                });

                // Collect results
                foreach (var collision in localCollisions)
                {
                    _collisions.Add(collision);
                }
            }
            else
            {
                foreach (var pair in _broadPhasePairs)
                {
                    var collision = CollisionDetector.Detect(pair.A, pair.B);
                    if (collision.HasCollision)
                    {
                        _collisions.Add(collision);
                    }
                }
            }
        }

        private void SolveConstraints(double dt)
        {
            _solver.VelocityIterations = SolverIterations;
            _solver.PositionIterations = PositionIterations;
            _solver.RestitutionThreshold = RestitutionThreshold;
            _solver.WarmStarting = UseWarmStarting;
            _solver.SolveVelocities(_collisions);
        }

        private void IntegratePositions(double dt)
        {
            bool useCcd = UseCCD;

            foreach (var body in _bodies)
            {
                if (body is not RigidBody rb || !rb.IsActive || rb.IsSleeping || rb.BodyType == BodyType.Static)
                    continue;

                double stepTime = dt;
                if (useCcd && rb.BodyType == BodyType.Dynamic && ContinuousCollision.NeedsSweep(rb, dt, out var swept))
                {
                    // Candidates come from the spatial hash built for this step; the
                    // margin covers bodies that have already been moved.
                    var candidates = _spatialHash.Query(swept.Expanded(0.5 * _spatialHash.CellSize));
                    stepTime = ContinuousCollision.TimeOfImpact(rb, candidates, dt, swept);
                }

                rb.IntegratePosition(stepTime);
            }

            // Remove the overlap that is left once the bodies have moved.
            _solver.SolvePositions();
        }

        #endregion

        #region CCD (Continuous Collision Detection)

        /// <summary>
        /// Sweeps a body along a displacement and reports the first body it would hit.
        /// The world uses the same test internally when <see cref="UseCCD"/> is enabled.
        /// </summary>
        /// <param name="body">The moving body.</param>
        /// <param name="displacement">World-space displacement to test.</param>
        /// <param name="firstHit">Information about the first hit, if any.</param>
        public bool RaycastCCD(IPhysicsBody body, Vector3D displacement, out CollisionInfo firstHit)
        {
            firstHit = default;

            double distance = displacement.Magnitude;
            if (distance < PhysicsConstants.Epsilon)
                return false;

            var direction = displacement / distance;
            double radius = 0;
            if (body is RigidBody rb)
            {
                radius = rb.ShapeType == CollisionShapeType.Sphere
                    ? rb.Radius
                    : Math.Min(rb.HalfExtents.X, Math.Min(rb.HalfExtents.Y, rb.HalfExtents.Z));
            }

            IPhysicsBody? hitBody = null;
            double hitTravel = distance;
            var start = body.Position;

            foreach (var other in _bodies)
            {
                if (ReferenceEquals(other, body) || !other.IsActive)
                    continue;

                // Conservative advancement of the body's core sphere.
                double travelled = 0;
                for (int iteration = 0; iteration < 32 && travelled < hitTravel; iteration++)
                {
                    double gap = CollisionDetector.DistanceToBody(start + direction * travelled, other) - radius;
                    if (gap <= 1e-4)
                    {
                        if (iteration > 0)
                        {
                            hitTravel = travelled;
                            hitBody = other;
                        }
                        break;
                    }
                    travelled += gap;
                }
            }

            if (hitBody == null)
                return false;

            firstHit = CollisionInfo.Create(body, hitBody, direction, 0, start + direction * (hitTravel + radius));
            return true;
        }

        #endregion

        #region Queries

        /// <summary>
        /// Casts a ray through the world using spatial hash acceleration.
        /// </summary>
        public RaycastHit? Raycast(Vector3D origin, Vector3D direction, double maxDistance = 1000)
        {
            direction = direction.Normalized;
            RaycastHit? closest = null;
            double closestDist = maxDistance;

            // Use spatial hash for faster raycast
            foreach (var body in _spatialHash.RaycastAll(origin, direction, maxDistance))
            {
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

            return closest;
        }

        /// <summary>
        /// Gets all bodies within a sphere using spatial hash.
        /// </summary>
        public List<IPhysicsBody> QuerySphere(Vector3D center, double radius)
        {
            var results = new List<IPhysicsBody>();
            double radiusSq = radius * radius;

            var queryAABB = new AABB(
                center - new Vector3D(radius, radius, radius),
                center + new Vector3D(radius, radius, radius)
            );

            foreach (var body in _spatialHash.Query(queryAABB))
            {
                if (Vector3D.DistanceSquared(center, body.Position) <= radiusSq)
                {
                    results.Add(body);
                }
            }

            return results;
        }

        /// <summary>
        /// Gets all bodies within an AABB using spatial hash.
        /// </summary>
        public List<IPhysicsBody> QueryAABB(AABB bounds)
        {
            var results = new List<IPhysicsBody>();

            foreach (var body in _spatialHash.Query(bounds))
            {
                if (bounds.Intersects(body.BoundingBox))
                {
                    results.Add(body);
                }
            }

            return results;
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

            foreach (var body in _bodies)
            {
                body.WakeUp();
            }
        }

        /// <summary>
        /// Gets performance statistics as a formatted string.
        /// </summary>
        public string GetPerformanceStats()
        {
            return $"Bodies: {_bodies.Count} (Dynamic: {_dynamicBodies.Count})\n" +
                   $"Islands: {_islandSolver.IslandCount} (Sleeping: {_islandSolver.SleepingIslandCount})\n" +
                   $"Broad Phase: {BroadPhaseTimeMs:F2}ms ({_broadPhasePairs.Count} pairs)\n" +
                   $"Narrow Phase: {NarrowPhaseTimeMs:F2}ms ({_collisions.Count} contacts)\n" +
                   $"Solver: {SolverTimeMs:F2}ms\n" +
                   $"Integration: {IntegrationTimeMs:F2}ms\n" +
                   $"Spatial Hash Cells: {_spatialHash.OccupiedCellCount}";
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

        #endregion
    }

}
