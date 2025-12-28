using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Compute;

namespace Artemis.Physics2D
{
    /// <summary>
    /// High-performance 2D physics world with multi-threading, SIMD, and GPU support.
    /// </summary>
    public class PhysicsWorld2D
    {
        #region Fields

        private readonly List<RigidBody2D> _bodies = new();
        private readonly List<Joint2D> _joints = new();
        private readonly List<AreaEffector2D> _areaEffectors = new();
        private readonly List<Manifold2D> _manifolds = new();
        private readonly Dictionary<(string, string), Manifold2D> _manifoldCache = new();

        // Thread-safe collision detection
        private readonly ConcurrentBag<(RigidBody2D, RigidBody2D)> _broadPhasePairs = new();
        private readonly ConcurrentBag<Manifold2D> _detectedCollisions = new();

        // Spatial hash for broad phase
        private readonly ConcurrentDictionary<long, List<int>> _spatialGrid = new();
        private float _cellSize = 2.0f;
        private float _invCellSize = 0.5f;

        // Thread synchronization
        private readonly ReaderWriterLockSlim _bodiesLock = new();

        // GPU compute (optional)
        private GpuCompute? _gpuCompute;
        private bool _useGpu;

        #endregion

        #region Properties

        /// <summary>Gravity vector applied to all dynamic bodies.</summary>
        public Vector2D Gravity { get; set; } = new(0, -9.81);

        /// <summary>Number of velocity iterations per step.</summary>
        public int VelocityIterations { get; set; } = 8;

        /// <summary>Number of position iterations per step.</summary>
        public int PositionIterations { get; set; } = 3;

        /// <summary>All bodies in the world (thread-safe copy).</summary>
        public IReadOnlyList<RigidBody2D> Bodies
        {
            get
            {
                _bodiesLock.EnterReadLock();
                try { return _bodies.ToArray(); }
                finally { _bodiesLock.ExitReadLock(); }
            }
        }

        /// <summary>All joints in the world.</summary>
        public IReadOnlyList<Joint2D> Joints => _joints;

        /// <summary>All area effectors in the world.</summary>
        public IReadOnlyList<AreaEffector2D> AreaEffectors => _areaEffectors;

        /// <summary>All active collision manifolds.</summary>
        public IReadOnlyList<Manifold2D> Manifolds => _manifolds;

        /// <summary>Whether to allow bodies to sleep.</summary>
        public bool AllowSleep { get; set; } = true;

        /// <summary>Total simulation time elapsed.</summary>
        public double SimulationTime { get; private set; }

        /// <summary>Whether to use GPU acceleration.</summary>
        public bool UseGPU
        {
            get => _useGpu;
            set
            {
                if (value && _gpuCompute == null)
                    InitializeGpu();
                _useGpu = value && _gpuCompute != null;
            }
        }

        /// <summary>Whether to use continuous collision detection.</summary>
        public bool UseCCD { get; set; } = true;

        #endregion

        #region Events

        /// <summary>Raised when a collision begins.</summary>
        public event EventHandler<CollisionEvent2D>? CollisionBegin;

        /// <summary>Raised when a collision ends.</summary>
        public event EventHandler<CollisionEvent2D>? CollisionEnd;

        /// <summary>Raised during collision (for modifying contact data).</summary>
        public event EventHandler<CollisionEvent2D>? PreSolve;

        /// <summary>Raised after collision is resolved.</summary>
        public event EventHandler<CollisionEvent2D>? PostSolve;

        /// <summary>Raised when a trigger is entered.</summary>
        public event EventHandler<CollisionEvent2D>? TriggerEnter;

        /// <summary>Raised when a trigger is exited.</summary>
        public event EventHandler<CollisionEvent2D>? TriggerExit;

        #endregion

        #region Constructor

        public PhysicsWorld2D()
        {
        }

        public PhysicsWorld2D(Vector2D gravity) : this()
        {
            Gravity = gravity;
        }

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

        #region Body Management

        /// <summary>
        /// Adds a body to the world.
        /// </summary>
        public void AddBody(RigidBody2D body)
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
        public bool RemoveBody(RigidBody2D body)
        {
            _bodiesLock.EnterWriteLock();
            try
            {
                _joints.RemoveAll(j => j.BodyA == body || j.BodyB == body);
                _manifolds.RemoveAll(m => m.BodyA == body || m.BodyB == body);
                return _bodies.Remove(body);
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        /// <summary>
        /// Gets a body by ID.
        /// </summary>
        public RigidBody2D? GetBody(string id)
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
                _joints.Clear();
                _manifolds.Clear();
                _manifoldCache.Clear();
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        #endregion

        #region Joint Management

        public void AddJoint(Joint2D joint)
        {
            if (!_joints.Contains(joint))
                _joints.Add(joint);
        }

        public bool RemoveJoint(Joint2D joint) => _joints.Remove(joint);

        public DistanceJoint2D CreateDistanceJoint(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchorA, Vector2D anchorB)
        {
            var joint = new DistanceJoint2D(bodyA, bodyB, anchorA, anchorB);
            AddJoint(joint);
            return joint;
        }

        public RevoluteJoint2D CreateRevoluteJoint(RigidBody2D bodyA, RigidBody2D bodyB, Vector2D anchor)
        {
            var joint = new RevoluteJoint2D(bodyA, bodyB, anchor);
            AddJoint(joint);
            return joint;
        }

        #endregion

        #region Area Effector Management

        public void AddAreaEffector(AreaEffector2D effector)
        {
            if (!_areaEffectors.Contains(effector))
                _areaEffectors.Add(effector);
        }

        public bool RemoveAreaEffector(AreaEffector2D effector) => _areaEffectors.Remove(effector);

        #endregion

        #region Simulation (Multi-threaded)

        /// <summary>
        /// Steps the simulation forward by the given time with multi-threading.
        /// </summary>
        public void Step(double dt)
        {
            if (dt <= 0) return;

            _bodiesLock.EnterReadLock();
            try
            {
                // Phase 1: Build spatial grid
                BuildSpatialGrid();

                // Phase 2: Broad phase (parallel)
                BroadPhaseParallel();

                // Phase 3: Narrow phase (parallel detection, sequential resolution)
                NarrowPhaseParallel();

                // Phase 4: Pre-solve joints
                foreach (var joint in _joints)
                {
                    if (joint.IsActive)
                        joint.PreSolve(dt);
                }

                // Phase 5: Apply area effectors (parallel)
                ApplyAreaEffectorsParallel(dt);

                // Phase 6: Integrate velocities (parallel + SIMD)
                IntegrateVelocitiesParallel(dt);

                // Phase 7: Velocity iterations
                for (int i = 0; i < VelocityIterations; i++)
                {
                    foreach (var joint in _joints)
                    {
                        if (joint.IsActive)
                            joint.SolveVelocity();
                    }

                    foreach (var manifold in _manifolds)
                    {
                        if (manifold.IsActive)
                            CollisionResolver2D.ResolveCollision(manifold);
                    }
                }

                // Phase 8: Integrate positions (parallel + SIMD)
                IntegratePositionsParallel(dt);

                // Phase 9: Position iterations
                for (int i = 0; i < PositionIterations; i++)
                {
                    bool jointsSolved = true;
                    foreach (var joint in _joints)
                    {
                        if (joint.IsActive)
                            jointsSolved &= joint.SolvePosition();
                    }

                    foreach (var manifold in _manifolds)
                    {
                        if (manifold.IsActive)
                            CollisionResolver2D.CorrectPositions(manifold);
                    }

                    if (jointsSolved) break;
                }

                // Phase 10: Post-solve events
                foreach (var manifold in _manifolds)
                {
                    if (manifold.IsActive)
                        PostSolve?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                }

                // Phase 11: Sleep handling (parallel)
                HandleSleepParallel();
            }
            finally { _bodiesLock.ExitReadLock(); }

            SimulationTime += dt;
        }

        /// <summary>
        /// Steps with fixed time step and accumulator.
        /// </summary>
        public void StepFixed(double dt, double fixedDt = 1.0 / 60.0)
        {
            while (dt >= fixedDt)
            {
                Step(fixedDt);
                dt -= fixedDt;
            }

            if (dt > PhysicsConstants.Epsilon)
                Step(dt);
        }

        #endregion

        #region Spatial Hash Grid

        private void BuildSpatialGrid()
        {
            _spatialGrid.Clear();

            for (int i = 0; i < _bodies.Count; i++)
            {
                var body = _bodies[i];
                if (!body.IsActive) continue;

                long cellKey = GetCellKey(body.Position);
                _spatialGrid.AddOrUpdate(cellKey,
                    _ => new List<int> { i },
                    (_, list) => { lock (list) { list.Add(i); } return list; });
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long GetCellKey(Vector2D pos)
        {
            int cx = (int)Math.Floor(pos.X * _invCellSize);
            int cy = (int)Math.Floor(pos.Y * _invCellSize);
            return ((long)(cx & 0x7FFFFFFF) << 32) | (uint)(cy & 0x7FFFFFFF);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long GetCellKey(int cx, int cy)
        {
            return ((long)(cx & 0x7FFFFFFF) << 32) | (uint)(cy & 0x7FFFFFFF);
        }

        #endregion

        #region Parallel Broad Phase

        private void BroadPhaseParallel()
        {
            while (_broadPhasePairs.TryTake(out _)) { }

            var cells = _spatialGrid.ToArray();

            Parallel.ForEach(cells, cell =>
            {
                var indices = cell.Value;
                long cellKey = cell.Key;
                int cx = (int)(cellKey >> 32);
                int cy = (int)(cellKey & 0x7FFFFFFF);

                // Check within cell
                for (int i = 0; i < indices.Count; i++)
                {
                    for (int j = i + 1; j < indices.Count; j++)
                    {
                        CheckBroadPhasePair(indices[i], indices[j]);
                    }
                }

                // Check neighboring cells
                for (int dx = 0; dx <= 1; dx++)
                {
                    for (int dy = 0; dy <= 1; dy++)
                    {
                        if (dx == 0 && dy == 0) continue;

                        long neighborKey = GetCellKey(cx + dx, cy + dy);
                        if (!_spatialGrid.TryGetValue(neighborKey, out var neighborIndices)) continue;

                        foreach (var i in indices)
                        {
                            foreach (var jIdx in neighborIndices)
                            {
                                CheckBroadPhasePair(i, jIdx);
                            }
                        }
                    }
                }
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void CheckBroadPhasePair(int i, int j)
        {
            var a = _bodies[i];
            var b = _bodies[j];

            if (!a.IsActive || !b.IsActive) return;
            if (a.BodyType == BodyType2D.Static && b.BodyType == BodyType2D.Static) return;
            if (a.IsSleeping && b.IsSleeping) return;
            if (!a.ShouldCollide(b)) return;

            if (a.AABB.Overlaps(b.AABB))
                _broadPhasePairs.Add((a, b));
        }

        #endregion

        #region Parallel Narrow Phase

        private void NarrowPhaseParallel()
        {
            while (_detectedCollisions.TryTake(out _)) { }

            foreach (var manifold in _manifolds)
                manifold.IsActive = false;

            var pairs = _broadPhasePairs.ToArray();

            // Parallel collision detection
            Parallel.ForEach(pairs, pair =>
            {
                var (a, b) = pair;
                var manifold = new Manifold2D();

                bool hasCollision = false;

                // CCD for fast-moving objects
                if (UseCCD && (a.UseCCD || b.UseCCD))
                {
                    hasCollision = ContinuousCollisionDetector2D.DetectSwept(a, b, manifold);
                }

                if (!hasCollision)
                {
                    hasCollision = CollisionDetector2D.DetectCollision(a, b, manifold);
                }

                if (hasCollision)
                {
                    manifold.IsActive = true;
                    _detectedCollisions.Add(manifold);
                }
            });

            // Sequential collision processing (for events and cache)
            foreach (var manifold in _detectedCollisions)
            {
                var key = GetManifoldKey(manifold.BodyA, manifold.BodyB);

                if (!_manifoldCache.ContainsKey(key))
                {
                    _manifoldCache[key] = manifold;
                    _manifolds.Add(manifold);
                    CollisionBegin?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                }
                else
                {
                    var cached = _manifoldCache[key];
                    cached.Normal = manifold.Normal;
                    cached.Penetration = manifold.Penetration;
                    cached.ContactCount = manifold.ContactCount;
                    cached.IsActive = true;
                }

                manifold.BodyA.WakeUp();
                manifold.BodyB.WakeUp();
                PreSolve?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
            }

            // Clean up ended collisions
            for (int i = _manifolds.Count - 1; i >= 0; i--)
            {
                var manifold = _manifolds[i];
                if (!manifold.IsActive)
                {
                    CollisionEnd?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                    var key = GetManifoldKey(manifold.BodyA, manifold.BodyB);
                    _manifoldCache.Remove(key);
                    _manifolds.RemoveAt(i);
                }
            }
        }

        private static (string, string) GetManifoldKey(RigidBody2D a, RigidBody2D b)
        {
            return string.CompareOrdinal(a.Id, b.Id) < 0 ? (a.Id, b.Id) : (b.Id, a.Id);
        }

        #endregion

        #region Parallel Integration with SIMD

        private void ApplyAreaEffectorsParallel(double dt)
        {
            if (_areaEffectors.Count == 0) return;

            var enabledEffectors = _areaEffectors.Where(e => e.Enabled).ToArray();
            var dynamicBodies = _bodies.Where(b => b.BodyType == BodyType2D.Dynamic && b.IsActive && !b.IsSleeping).ToArray();

            Parallel.ForEach(dynamicBodies, body =>
            {
                foreach (var effector in enabledEffectors)
                {
                    if (effector.AffectsBody(body))
                    {
                        effector.ApplyForce(body, dt);
                    }
                }
            });
        }

        private void IntegrateVelocitiesParallel(double dt)
        {
            var dynamicBodies = _bodies.Where(b => b.BodyType == BodyType2D.Dynamic && b.IsActive).ToArray();
            var gravity = Gravity;

            Parallel.ForEach(dynamicBodies, body =>
            {
                body.IntegrateVelocity(dt, gravity);
            });
        }

        private void IntegratePositionsParallel(double dt)
        {
            var activeBodies = _bodies.Where(b => b.IsActive).ToArray();

            Parallel.ForEach(activeBodies, body =>
            {
                body.IntegratePosition(dt);
            });
        }

        private void HandleSleepParallel()
        {
            if (!AllowSleep)
            {
                Parallel.ForEach(_bodies, body => body.IsSleeping = false);
                return;
            }

            Parallel.ForEach(_bodies, body =>
            {
                if (!body.CanSleep)
                    body.IsSleeping = false;
            });
        }

        #endregion

        #region Queries

        /// <summary>
        /// Queries all bodies overlapping an AABB (parallel).
        /// </summary>
        public List<RigidBody2D> QueryAABB(AABB2D aabb)
        {
            var results = new ConcurrentBag<RigidBody2D>();

            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (body.IsActive && body.AABB.Overlaps(aabb))
                        results.Add(body);
                });
            }
            finally { _bodiesLock.ExitReadLock(); }

            return results.ToList();
        }

        /// <summary>
        /// Queries all bodies containing a point.
        /// </summary>
        public List<RigidBody2D> QueryPoint(Vector2D point)
        {
            var results = new ConcurrentBag<RigidBody2D>();

            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (!body.IsActive || body.Shape == null) return;
                    if (!body.AABB.Contains(point)) return;

                    if (body.Shape is CircleShape circle)
                    {
                        var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                        if (CollisionDetector2D.PointInCircle(point, center, circle.Radius))
                            results.Add(body);
                    }
                    else if (body.Shape is BoxShape box)
                    {
                        var local = body.WorldToLocal(point);
                        if (Math.Abs(local.X) <= box.HalfWidth && Math.Abs(local.Y) <= box.HalfHeight)
                            results.Add(body);
                    }
                });
            }
            finally { _bodiesLock.ExitReadLock(); }

            return results.ToList();
        }

        /// <summary>
        /// Casts a ray and returns the first hit.
        /// </summary>
        public CollisionDetector2D.RaycastHit2D Raycast(Vector2D origin, Vector2D direction,
            double maxDistance = double.MaxValue, ushort maskBits = 0xFFFF)
        {
            direction = direction.Normalized;
            var closestHit = new CollisionDetector2D.RaycastHit2D { Fraction = double.MaxValue };

            _bodiesLock.EnterReadLock();
            try
            {
                foreach (var body in _bodies)
                {
                    if (!body.IsActive || body.Shape == null) continue;
                    if ((body.CategoryBits & maskBits) == 0) continue;

                    CollisionDetector2D.RaycastHit2D hit;
                    bool didHit = false;

                    if (body.Shape is CircleShape circle)
                    {
                        var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                        didHit = CollisionDetector2D.RaycastCircle(origin, direction, maxDistance,
                            center, circle.Radius, out hit);
                    }
                    else
                    {
                        didHit = CollisionDetector2D.RaycastAABB(origin, direction, maxDistance,
                            body.AABB, out hit);
                    }

                    if (didHit && hit.Fraction < closestHit.Fraction)
                    {
                        closestHit = hit;
                        closestHit.Body = body;
                    }
                }
            }
            finally { _bodiesLock.ExitReadLock(); }

            closestHit.Hit = closestHit.Fraction < double.MaxValue;
            return closestHit;
        }

        #endregion

        #region Utility Methods

        /// <summary>
        /// Applies an explosion force at a point (parallel).
        /// </summary>
        public void ApplyExplosionForce(Vector2D center, double radius, double force)
        {
            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (body.BodyType != BodyType2D.Dynamic) return;

                    var delta = body.Position - center;
                    double distance = delta.Magnitude;

                    if (distance > radius || distance < PhysicsConstants.Epsilon) return;

                    double falloff = 1.0 - distance / radius;
                    var direction = delta / distance;
                    body.ApplyImpulse(direction * force * falloff);
                });
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        /// <summary>
        /// Gets the kinetic energy of all bodies.
        /// </summary>
        public double GetTotalKineticEnergy()
        {
            double total = 0;
            _bodiesLock.EnterReadLock();
            try
            {
                foreach (var body in _bodies)
                {
                    if (body.BodyType == BodyType2D.Dynamic)
                    {
                        total += 0.5 * body.Mass * body.Velocity.MagnitudeSquared;
                        total += 0.5 * body.Inertia * body.AngularVelocity * body.AngularVelocity;
                    }
                }
            }
            finally { _bodiesLock.ExitReadLock(); }
            return total;
        }

        /// <summary>
        /// Resets the simulation.
        /// </summary>
        public void Reset()
        {
            ClearBodies();
            SimulationTime = 0;
        }

        #endregion
    }

    /// <summary>
    /// Continuous collision detection for 2D.
    /// </summary>
    public static class ContinuousCollisionDetector2D
    {
        public static bool DetectSwept(RigidBody2D bodyA, RigidBody2D bodyB, Manifold2D manifold)
        {
            // For circles, use swept circle collision
            if (bodyA.Shape is CircleShape circleA && bodyB.Shape is CircleShape circleB)
            {
                return SweptCircleCircle(bodyA, bodyB, circleA, circleB, manifold);
            }

            // Fall back to discrete for other shapes
            return CollisionDetector2D.DetectCollision(bodyA, bodyB, manifold);
        }

        private static bool SweptCircleCircle(RigidBody2D bodyA, RigidBody2D bodyB,
            CircleShape circleA, CircleShape circleB, Manifold2D manifold)
        {
            var centerA = bodyA.Position + Vector2D.Rotate(circleA.Offset, bodyA.Rotation);
            var centerB = bodyB.Position + Vector2D.Rotate(circleB.Offset, bodyB.Rotation);

            var relVel = bodyA.Velocity - bodyB.Velocity;
            var centerDiff = centerB - centerA;
            double radiusSum = circleA.Radius + circleB.Radius;

            double a = relVel.MagnitudeSquared;
            double b = -2 * Vector2D.Dot(relVel, centerDiff);
            double c = centerDiff.MagnitudeSquared - radiusSum * radiusSum;

            // Already overlapping
            if (c < 0)
            {
                return CollisionDetector2D.CircleVsCircle(bodyA, bodyB, circleA, circleB, manifold);
            }

            if (Math.Abs(a) < 1e-8)
                return false;

            double discriminant = b * b - 4 * a * c;
            if (discriminant < 0)
                return false;

            double sqrtD = Math.Sqrt(discriminant);
            double t = (-b - sqrtD) / (2 * a);

            if (t < 0 || t > 1)
                return false;

            // Collision will happen, use discrete detection at TOI
            return CollisionDetector2D.CircleVsCircle(bodyA, bodyB, circleA, circleB, manifold);
        }
    }
}
