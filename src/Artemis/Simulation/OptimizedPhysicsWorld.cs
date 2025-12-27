using System;
using System.Collections.Generic;
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
        private const double SleepAngularThreshold = 0.05;
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

        // Warm starting cache (stores impulses from previous frame)
        private readonly Dictionary<(IPhysicsBody, IPhysicsBody), WarmStartData> _warmStartCache;

        // Threading
        private readonly int _threadCount;
        private readonly ParallelOptions _parallelOptions;

        // Timing
        private double _accumulator;
        private readonly Stopwatch _profileTimer;

        #endregion

        #region Nested Types

        private struct WarmStartData
        {
            public Vector3D NormalImpulse;
            public Vector3D TangentImpulse;
            public Vector3D ContactPoint;
        }

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
            _warmStartCache = new Dictionary<(IPhysicsBody, IPhysicsBody), WarmStartData>(256);

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
            return _bodies.Remove(body);
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
            _warmStartCache.Clear();
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

            // 5. Build islands
            _islandSolver.BuildIslands(_bodies, _collisions);

            // 6. Solve constraints
            _profileTimer.Restart();
            SolveConstraints(dt);
            SolverTimeMs = _profileTimer.Elapsed.TotalMilliseconds;

            // 7. Integrate positions
            IntegratePositions(dt);

            // 8. Update sleeping
            UpdateSleeping(dt);

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
            if (UseMultiThreading && _dynamicBodies.Count > 100)
            {
                Parallel.ForEach(_dynamicBodies, _parallelOptions, body =>
                {
                    ApplyForcesToBody(body);
                });
            }
            else
            {
                foreach (var body in _dynamicBodies)
                {
                    ApplyForcesToBody(body);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void ApplyForcesToBody(IPhysicsBody body)
        {
            if (!body.IsActive || body.IsSleeping)
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
            if (UseMultiThreading && _dynamicBodies.Count > 100)
            {
                Parallel.ForEach(_dynamicBodies, _parallelOptions, body =>
                {
                    if (body.IsActive && !body.IsSleeping)
                    {
                        IntegrateVelocity(body, dt);
                    }
                });
            }
            else
            {
                foreach (var body in _dynamicBodies)
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
            if (body is RigidBody rb)
            {
                // Semi-implicit Euler for velocity
                rb.Velocity += rb.Acceleration * dt;
                rb.AngularVelocity += rb.AngularAcceleration * dt;

                // Apply damping
                rb.Velocity *= Math.Pow(1.0 - rb.LinearDamping, dt);
                rb.AngularVelocity *= Math.Pow(1.0 - rb.AngularDamping, dt);
            }
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
            if (UseMultiThreading)
            {
                // Solve islands in parallel
                _islandSolver.SolveParallel(island =>
                {
                    SolveIsland(island, dt);
                }, _threadCount);
            }
            else
            {
                _islandSolver.SolveSequential(island =>
                {
                    SolveIsland(island, dt);
                });
            }
        }

        private void SolveIsland(IslandSolver.Island island, double dt)
        {
            // Velocity solver iterations
            for (int i = 0; i < SolverIterations; i++)
            {
                foreach (var collision in island.Collisions)
                {
                    SolveVelocityConstraint(collision, dt);
                }
            }

            // Position solver iterations
            for (int i = 0; i < PositionIterations; i++)
            {
                foreach (var collision in island.Collisions)
                {
                    SolvePositionConstraint(collision);
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void SolveVelocityConstraint(CollisionInfo collision, double dt)
        {
            var bodyA = collision.BodyA as RigidBody;
            var bodyB = collision.BodyB as RigidBody;
            if (bodyA == null && bodyB == null) return;

            // Get masses
            double invMassA = bodyA?.BodyType == BodyType.Dynamic ? 1.0 / bodyA.Mass : 0;
            double invMassB = bodyB?.BodyType == BodyType.Dynamic ? 1.0 / bodyB.Mass : 0;

            if (invMassA == 0 && invMassB == 0) return;

            // Contact point relative to centers
            var rA = collision.ContactPoint - (bodyA?.Position ?? collision.ContactPoint);
            var rB = collision.ContactPoint - (bodyB?.Position ?? collision.ContactPoint);

            // Relative velocity at contact
            var vA = (bodyA?.Velocity ?? Vector3D.Zero) +
                     Vector3D.Cross(bodyA?.AngularVelocity ?? Vector3D.Zero, rA);
            var vB = (bodyB?.Velocity ?? Vector3D.Zero) +
                     Vector3D.Cross(bodyB?.AngularVelocity ?? Vector3D.Zero, rB);
            var relVel = vB - vA;

            // Normal velocity
            double normalVel = Vector3D.Dot(relVel, collision.Normal);

            // Don't resolve if separating
            if (normalVel > 0) return;

            // Restitution
            double restitution = Math.Min(
                bodyA?.Material.Restitution ?? 0.5,
                bodyB?.Material.Restitution ?? 0.5
            );

            // Only apply restitution above threshold
            if (Math.Abs(normalVel) < RestitutionThreshold)
            {
                restitution = 0;
            }

            // Effective mass
            var rAxN = Vector3D.Cross(rA, collision.Normal);
            var rBxN = Vector3D.Cross(rB, collision.Normal);

            double kNormal = invMassA + invMassB;
            if (bodyA != null)
            {
                kNormal += Vector3D.Dot(
                    Vector3D.Cross(bodyA.InverseInertiaTensorWorld * rAxN, rA),
                    collision.Normal);
            }
            if (bodyB != null)
            {
                kNormal += Vector3D.Dot(
                    Vector3D.Cross(bodyB.InverseInertiaTensorWorld * rBxN, rB),
                    collision.Normal);
            }

            // Normal impulse
            double j = -(1 + restitution) * normalVel / kNormal;

            // Warm starting
            var pairKey = (collision.BodyA, collision.BodyB);
            if (UseWarmStarting && _warmStartCache.TryGetValue(pairKey, out var warmStart))
            {
                // Apply cached impulse scaled down
                var warmImpulse = warmStart.NormalImpulse * 0.8;
                if (bodyA != null && bodyA.BodyType == BodyType.Dynamic)
                {
                    bodyA.Velocity -= warmImpulse * invMassA;
                    bodyA.AngularVelocity -= bodyA.InverseInertiaTensorWorld * Vector3D.Cross(rA, warmImpulse);
                }
                if (bodyB != null && bodyB.BodyType == BodyType.Dynamic)
                {
                    bodyB.Velocity += warmImpulse * invMassB;
                    bodyB.AngularVelocity += bodyB.InverseInertiaTensorWorld * Vector3D.Cross(rB, warmImpulse);
                }
            }

            // Apply normal impulse
            var impulse = collision.Normal * j;

            if (bodyA != null && bodyA.BodyType == BodyType.Dynamic)
            {
                bodyA.Velocity -= impulse * invMassA;
                bodyA.AngularVelocity -= bodyA.InverseInertiaTensorWorld * Vector3D.Cross(rA, impulse);
            }
            if (bodyB != null && bodyB.BodyType == BodyType.Dynamic)
            {
                bodyB.Velocity += impulse * invMassB;
                bodyB.AngularVelocity += bodyB.InverseInertiaTensorWorld * Vector3D.Cross(rB, impulse);
            }

            // Friction
            var tangent = relVel - collision.Normal * normalVel;
            double tangentLen = tangent.Length;
            if (tangentLen > 1e-6)
            {
                tangent /= tangentLen;

                double friction = Math.Sqrt(
                    (bodyA?.Material.Friction ?? 0.5) *
                    (bodyB?.Material.Friction ?? 0.5)
                );

                double jt = -Vector3D.Dot(relVel, tangent) / kNormal;
                jt = Math.Clamp(jt, -Math.Abs(j) * friction, Math.Abs(j) * friction);

                var frictionImpulse = tangent * jt;

                if (bodyA != null && bodyA.BodyType == BodyType.Dynamic)
                {
                    bodyA.Velocity -= frictionImpulse * invMassA;
                }
                if (bodyB != null && bodyB.BodyType == BodyType.Dynamic)
                {
                    bodyB.Velocity += frictionImpulse * invMassB;
                }
            }

            // Store for warm starting
            _warmStartCache[pairKey] = new WarmStartData
            {
                NormalImpulse = impulse,
                ContactPoint = collision.ContactPoint
            };
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void SolvePositionConstraint(CollisionInfo collision)
        {
            if (collision.PenetrationDepth < 0.001)
                return;

            var bodyA = collision.BodyA as RigidBody;
            var bodyB = collision.BodyB as RigidBody;

            double invMassA = bodyA?.BodyType == BodyType.Dynamic ? 1.0 / bodyA.Mass : 0;
            double invMassB = bodyB?.BodyType == BodyType.Dynamic ? 1.0 / bodyB.Mass : 0;

            if (invMassA == 0 && invMassB == 0) return;

            // Baumgarte stabilization
            const double baumgarte = 0.2;
            const double slop = 0.005;

            double correction = baumgarte * Math.Max(collision.PenetrationDepth - slop, 0) /
                               (invMassA + invMassB);

            var correctionVec = collision.Normal * correction;

            if (bodyA != null && bodyA.BodyType == BodyType.Dynamic)
            {
                bodyA.Position -= correctionVec * invMassA;
            }
            if (bodyB != null && bodyB.BodyType == BodyType.Dynamic)
            {
                bodyB.Position += correctionVec * invMassB;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void IntegratePositions(double dt)
        {
            foreach (var body in _dynamicBodies)
            {
                if (!body.IsActive || body.IsSleeping)
                    continue;

                if (body is RigidBody rb)
                {
                    rb.Position += rb.Velocity * dt;

                    // Quaternion integration for rotation
                    var w = new Quaternion(0, rb.AngularVelocity.X, rb.AngularVelocity.Y, rb.AngularVelocity.Z);
                    var qDot = 0.5 * w * rb.Rotation;
                    rb.Rotation = (rb.Rotation + qDot * dt).Normalized();

                    // Clear forces
                    rb.ClearForces();
                }
            }
        }

        private void UpdateSleeping(double dt)
        {
            foreach (var body in _dynamicBodies)
            {
                if (!body.IsActive) continue;

                if (body is RigidBody rb)
                {
                    double linearSpeed = rb.Velocity.Length;
                    double angularSpeed = rb.AngularVelocity.Length;

                    if (linearSpeed < SleepVelocityThreshold && angularSpeed < SleepAngularThreshold)
                    {
                        rb.SleepTimer += dt;
                        if (rb.SleepTimer >= SleepTimeRequired)
                        {
                            rb.Sleep();
                        }
                    }
                    else
                    {
                        rb.SleepTimer = 0;
                        if (rb.IsSleeping)
                        {
                            rb.WakeUp();
                        }
                    }
                }
            }
        }

        #endregion

        #region CCD (Continuous Collision Detection)

        /// <summary>
        /// Performs continuous collision detection for fast-moving objects.
        /// </summary>
        public bool RaycastCCD(IPhysicsBody body, Vector3D displacement, out CollisionInfo firstHit)
        {
            firstHit = default;
            if (!UseCCD) return false;

            double speed = displacement.Length;
            if (speed < body.BoundingBox.Size.Length * 0.5)
            {
                return false; // Not moving fast enough
            }

            var direction = displacement / speed;

            // Expand AABB along motion path
            var sweepAABB = body.BoundingBox;
            var endAABB = new AABB(
                body.BoundingBox.Min + displacement,
                body.BoundingBox.Max + displacement
            );
            sweepAABB = AABB.Union(sweepAABB, endAABB);

            // Query spatial hash for potential hits
            var candidates = _spatialHash.Query(sweepAABB);

            double minT = double.MaxValue;

            foreach (var other in candidates)
            {
                if (other == body || !other.IsActive) continue;

                // Simple sphere-cast for CCD
                if (RaycastBody(body.Position, direction, speed, other, out double t, out var hitInfo))
                {
                    if (t < minT)
                    {
                        minT = t;
                        firstHit = hitInfo;
                    }
                }
            }

            return minT < double.MaxValue;
        }

        private bool RaycastBody(Vector3D origin, Vector3D direction, double maxDist,
            IPhysicsBody target, out double t, out CollisionInfo collision)
        {
            t = double.MaxValue;
            collision = default;

            // Simplified: just use AABB for now
            var aabb = target.BoundingBox;

            double tmin = 0;
            double tmax = maxDist;

            for (int i = 0; i < 3; i++)
            {
                double o = i == 0 ? origin.X : (i == 1 ? origin.Y : origin.Z);
                double d = i == 0 ? direction.X : (i == 1 ? direction.Y : direction.Z);
                double min = i == 0 ? aabb.Min.X : (i == 1 ? aabb.Min.Y : aabb.Min.Z);
                double max = i == 0 ? aabb.Max.X : (i == 1 ? aabb.Max.Y : aabb.Max.Z);

                if (Math.Abs(d) < 1e-8)
                {
                    if (o < min || o > max) return false;
                }
                else
                {
                    double t1 = (min - o) / d;
                    double t2 = (max - o) / d;
                    if (t1 > t2) (t1, t2) = (t2, t1);
                    tmin = Math.Max(tmin, t1);
                    tmax = Math.Min(tmax, t2);
                    if (tmin > tmax) return false;
                }
            }

            t = tmin;
            return true;
        }

        #endregion

        #region Queries

        /// <summary>
        /// Casts a ray through the world using spatial hash acceleration.
        /// </summary>
        public RaycastHit? Raycast(Vector3D origin, Vector3D direction, double maxDistance = 1000)
        {
            direction = direction.Normalized();
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
            _warmStartCache.Clear();

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

    // Thread-safe bag for parallel collection
    internal class ConcurrentBag<T>
    {
        private readonly List<T>[] _buckets;
        private readonly object[] _locks;

        public ConcurrentBag()
        {
            int count = Environment.ProcessorCount;
            _buckets = new List<T>[count];
            _locks = new object[count];
            for (int i = 0; i < count; i++)
            {
                _buckets[i] = new List<T>();
                _locks[i] = new object();
            }
        }

        public void Add(T item)
        {
            int bucket = Thread.CurrentThread.ManagedThreadId % _buckets.Length;
            lock (_locks[bucket])
            {
                _buckets[bucket].Add(item);
            }
        }

        public IEnumerator<T> GetEnumerator()
        {
            foreach (var bucket in _buckets)
            {
                foreach (var item in bucket)
                {
                    yield return item;
                }
            }
        }
    }
}
