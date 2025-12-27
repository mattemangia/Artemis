using System;
using System.Collections.Generic;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Event arguments for collision events.
    /// </summary>
    public class CollisionEvent2D : EventArgs
    {
        public RigidBody2D BodyA { get; }
        public RigidBody2D BodyB { get; }
        public Manifold2D Manifold { get; }

        public CollisionEvent2D(RigidBody2D bodyA, RigidBody2D bodyB, Manifold2D manifold)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            Manifold = manifold;
        }
    }

    /// <summary>
    /// 2D physics simulation world.
    /// </summary>
    public class PhysicsWorld2D
    {
        #region Fields

        private readonly List<RigidBody2D> _bodies = new();
        private readonly List<Joint2D> _joints = new();
        private readonly List<Manifold2D> _manifolds = new();
        private readonly Dictionary<(string, string), Manifold2D> _manifoldCache = new();
        private readonly List<(RigidBody2D, RigidBody2D)> _broadPhasePairs = new();

        #endregion

        #region Properties

        /// <summary>Gravity vector applied to all dynamic bodies.</summary>
        public Vector2D Gravity { get; set; } = new(0, -9.81);

        /// <summary>Number of velocity iterations per step.</summary>
        public int VelocityIterations { get; set; } = 8;

        /// <summary>Number of position iterations per step.</summary>
        public int PositionIterations { get; set; } = 3;

        /// <summary>All bodies in the world.</summary>
        public IReadOnlyList<RigidBody2D> Bodies => _bodies;

        /// <summary>All joints in the world.</summary>
        public IReadOnlyList<Joint2D> Joints => _joints;

        /// <summary>All active collision manifolds.</summary>
        public IReadOnlyList<Manifold2D> Manifolds => _manifolds;

        /// <summary>Whether to allow bodies to sleep.</summary>
        public bool AllowSleep { get; set; } = true;

        /// <summary>Total simulation time elapsed.</summary>
        public double SimulationTime { get; private set; }

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

        #endregion

        #region Body Management

        /// <summary>
        /// Adds a body to the world.
        /// </summary>
        public void AddBody(RigidBody2D body)
        {
            if (!_bodies.Contains(body))
            {
                _bodies.Add(body);
            }
        }

        /// <summary>
        /// Removes a body from the world.
        /// </summary>
        public bool RemoveBody(RigidBody2D body)
        {
            // Remove any joints connected to this body
            _joints.RemoveAll(j => j.BodyA == body || j.BodyB == body);

            // Remove manifolds involving this body
            _manifolds.RemoveAll(m => m.BodyA == body || m.BodyB == body);

            return _bodies.Remove(body);
        }

        /// <summary>
        /// Gets a body by ID.
        /// </summary>
        public RigidBody2D? GetBody(string id)
            => _bodies.Find(b => b.Id == id);

        /// <summary>
        /// Clears all bodies from the world.
        /// </summary>
        public void ClearBodies()
        {
            _bodies.Clear();
            _joints.Clear();
            _manifolds.Clear();
            _manifoldCache.Clear();
        }

        #endregion

        #region Joint Management

        /// <summary>
        /// Adds a joint to the world.
        /// </summary>
        public void AddJoint(Joint2D joint)
        {
            if (!_joints.Contains(joint))
            {
                _joints.Add(joint);
            }
        }

        /// <summary>
        /// Removes a joint from the world.
        /// </summary>
        public bool RemoveJoint(Joint2D joint)
            => _joints.Remove(joint);

        /// <summary>
        /// Creates and adds a distance joint.
        /// </summary>
        public DistanceJoint2D CreateDistanceJoint(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchorA, Vector2D anchorB)
        {
            var joint = new DistanceJoint2D(bodyA, bodyB, anchorA, anchorB);
            AddJoint(joint);
            return joint;
        }

        /// <summary>
        /// Creates and adds a revolute joint.
        /// </summary>
        public RevoluteJoint2D CreateRevoluteJoint(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchor)
        {
            var joint = new RevoluteJoint2D(bodyA, bodyB, anchor);
            AddJoint(joint);
            return joint;
        }

        /// <summary>
        /// Creates and adds a weld joint.
        /// </summary>
        public WeldJoint2D CreateWeldJoint(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchor)
        {
            var joint = new WeldJoint2D(bodyA, bodyB, anchor);
            AddJoint(joint);
            return joint;
        }

        /// <summary>
        /// Creates and adds a mouse joint for dragging.
        /// </summary>
        public MouseJoint2D CreateMouseJoint(RigidBody2D body, Vector2D target)
        {
            var joint = new MouseJoint2D(body, target);
            AddJoint(joint);
            return joint;
        }

        /// <summary>
        /// Creates and adds a rope joint.
        /// </summary>
        public RopeJoint2D CreateRopeJoint(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchorA, Vector2D anchorB, double maxLength)
        {
            var joint = new RopeJoint2D(bodyA, bodyB, anchorA, anchorB, maxLength);
            AddJoint(joint);
            return joint;
        }

        #endregion

        #region Simulation

        /// <summary>
        /// Steps the simulation forward by the given time.
        /// </summary>
        public void Step(double dt)
        {
            if (dt <= 0)
                return;

            // Broad phase collision detection
            BroadPhase();

            // Narrow phase collision detection
            NarrowPhase();

            // Pre-solve joints
            foreach (var joint in _joints)
            {
                if (joint.IsActive)
                    joint.PreSolve(dt);
            }

            // Integrate velocities
            foreach (var body in _bodies)
            {
                if (body.IsActive && body.BodyType == BodyType2D.Dynamic)
                {
                    body.IntegrateVelocity(dt, Gravity);
                }
            }

            // Velocity iterations
            for (int i = 0; i < VelocityIterations; i++)
            {
                // Solve joints
                foreach (var joint in _joints)
                {
                    if (joint.IsActive)
                        joint.SolveVelocity();
                }

                // Solve collisions
                foreach (var manifold in _manifolds)
                {
                    if (manifold.IsActive)
                        CollisionResolver2D.ResolveCollision(manifold);
                }
            }

            // Integrate positions
            foreach (var body in _bodies)
            {
                if (body.IsActive)
                {
                    body.IntegratePosition(dt);
                }
            }

            // Position iterations
            for (int i = 0; i < PositionIterations; i++)
            {
                // Solve joint positions
                bool jointsSolved = true;
                foreach (var joint in _joints)
                {
                    if (joint.IsActive)
                        jointsSolved &= joint.SolvePosition();
                }

                // Correct penetration
                foreach (var manifold in _manifolds)
                {
                    if (manifold.IsActive)
                        CollisionResolver2D.CorrectPositions(manifold);
                }

                if (jointsSolved)
                    break;
            }

            // Post-solve events
            foreach (var manifold in _manifolds)
            {
                if (manifold.IsActive)
                {
                    PostSolve?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                }
            }

            // Handle sleep
            if (AllowSleep)
            {
                foreach (var body in _bodies)
                {
                    if (!body.CanSleep)
                        body.IsSleeping = false;
                }
            }
            else
            {
                foreach (var body in _bodies)
                {
                    body.IsSleeping = false;
                }
            }

            SimulationTime += dt;
        }

        /// <summary>
        /// Steps with fixed time step and accumulator.
        /// </summary>
        public void StepFixed(double dt, double fixedDt = 1.0 / 60.0)
        {
            // Simple approach: step as many times as needed
            while (dt >= fixedDt)
            {
                Step(fixedDt);
                dt -= fixedDt;
            }

            if (dt > PhysicsConstants.Epsilon)
            {
                Step(dt);
            }
        }

        #endregion

        #region Collision Detection

        private void BroadPhase()
        {
            _broadPhasePairs.Clear();

            // Simple O(n^2) broad phase with AABB checks
            for (int i = 0; i < _bodies.Count; i++)
            {
                var a = _bodies[i];
                if (!a.IsActive || a.Shape == null)
                    continue;

                for (int j = i + 1; j < _bodies.Count; j++)
                {
                    var b = _bodies[j];
                    if (!b.IsActive || b.Shape == null)
                        continue;

                    // Skip static-static pairs
                    if (a.BodyType == BodyType2D.Static && b.BodyType == BodyType2D.Static)
                        continue;

                    // Skip if both sleeping
                    if (a.IsSleeping && b.IsSleeping)
                        continue;

                    // Check collision filters
                    if (!a.ShouldCollide(b))
                        continue;

                    // Check connected by joint
                    if (!CanCollideThroughJoints(a, b))
                        continue;

                    // AABB overlap test
                    if (a.AABB.Overlaps(b.AABB))
                    {
                        _broadPhasePairs.Add((a, b));
                    }
                }
            }
        }

        private bool CanCollideThroughJoints(RigidBody2D a, RigidBody2D b)
        {
            foreach (var joint in _joints)
            {
                if (!joint.CollideConnected)
                {
                    if ((joint.BodyA == a && joint.BodyB == b) ||
                        (joint.BodyA == b && joint.BodyB == a))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private void NarrowPhase()
        {
            // Mark existing manifolds as potentially inactive
            foreach (var manifold in _manifolds)
            {
                manifold.IsActive = false;
            }

            foreach (var (a, b) in _broadPhasePairs)
            {
                var key = GetManifoldKey(a, b);

                // Try to get existing manifold
                if (!_manifoldCache.TryGetValue(key, out var manifold))
                {
                    manifold = new Manifold2D();
                    _manifoldCache[key] = manifold;
                    _manifolds.Add(manifold);
                }

                // Perform narrow phase
                if (CollisionDetector2D.DetectCollision(a, b, manifold))
                {
                    bool wasActive = manifold.IsActive;
                    manifold.IsActive = true;

                    // Wake up bodies
                    a.WakeUp();
                    b.WakeUp();

                    // Fire events
                    if (!wasActive)
                    {
                        CollisionBegin?.Invoke(this, new CollisionEvent2D(a, b, manifold));
                    }

                    PreSolve?.Invoke(this, new CollisionEvent2D(a, b, manifold));
                }
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

        #region Queries

        /// <summary>
        /// Queries all bodies overlapping an AABB.
        /// </summary>
        public List<RigidBody2D> QueryAABB(AABB2D aabb)
        {
            var results = new List<RigidBody2D>();
            foreach (var body in _bodies)
            {
                if (body.IsActive && body.AABB.Overlaps(aabb))
                {
                    results.Add(body);
                }
            }
            return results;
        }

        /// <summary>
        /// Queries all bodies containing a point.
        /// </summary>
        public List<RigidBody2D> QueryPoint(Vector2D point)
        {
            var results = new List<RigidBody2D>();
            foreach (var body in _bodies)
            {
                if (!body.IsActive || body.Shape == null)
                    continue;

                if (!body.AABB.Contains(point))
                    continue;

                // More precise check based on shape
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
                else if (body.Shape is PolygonShape poly)
                {
                    var worldVerts = new Vector2D[poly.VertexCount];
                    for (int i = 0; i < poly.VertexCount; i++)
                        worldVerts[i] = body.LocalToWorld(poly.Vertices[i]);

                    if (CollisionDetector2D.PointInPolygon(point, worldVerts))
                        results.Add(body);
                }
            }
            return results;
        }

        /// <summary>
        /// Casts a ray and returns the first hit.
        /// </summary>
        public CollisionDetector2D.RaycastHit2D Raycast(Vector2D origin, Vector2D direction,
            double maxDistance = double.MaxValue, ushort maskBits = 0xFFFF)
        {
            direction = direction.Normalized;
            var closestHit = new CollisionDetector2D.RaycastHit2D();
            closestHit.Fraction = double.MaxValue;

            foreach (var body in _bodies)
            {
                if (!body.IsActive || body.Shape == null)
                    continue;

                if ((body.CategoryBits & maskBits) == 0)
                    continue;

                // Quick AABB check
                var rayAABB = new AABB2D(
                    Vector2D.Min(origin, origin + direction * maxDistance),
                    Vector2D.Max(origin, origin + direction * maxDistance)
                );
                if (!body.AABB.Overlaps(rayAABB))
                    continue;

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

            closestHit.Hit = closestHit.Fraction < double.MaxValue;
            return closestHit;
        }

        /// <summary>
        /// Casts a ray and returns all hits.
        /// </summary>
        public List<CollisionDetector2D.RaycastHit2D> RaycastAll(Vector2D origin, Vector2D direction,
            double maxDistance = double.MaxValue, ushort maskBits = 0xFFFF)
        {
            var results = new List<CollisionDetector2D.RaycastHit2D>();
            direction = direction.Normalized;

            foreach (var body in _bodies)
            {
                if (!body.IsActive || body.Shape == null)
                    continue;

                if ((body.CategoryBits & maskBits) == 0)
                    continue;

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

                if (didHit)
                {
                    hit.Body = body;
                    results.Add(hit);
                }
            }

            results.Sort((a, b) => a.Fraction.CompareTo(b.Fraction));
            return results;
        }

        #endregion

        #region Utility Methods

        /// <summary>
        /// Applies an explosion force at a point.
        /// </summary>
        public void ApplyExplosionForce(Vector2D center, double radius, double force)
        {
            foreach (var body in _bodies)
            {
                if (body.BodyType != BodyType2D.Dynamic)
                    continue;

                var delta = body.Position - center;
                double distance = delta.Magnitude;

                if (distance > radius || distance < PhysicsConstants.Epsilon)
                    continue;

                double falloff = 1.0 - distance / radius;
                var direction = delta / distance;
                body.ApplyImpulse(direction * force * falloff);
            }
        }

        /// <summary>
        /// Gets the kinetic energy of all bodies.
        /// </summary>
        public double GetTotalKineticEnergy()
        {
            double total = 0;
            foreach (var body in _bodies)
            {
                if (body.BodyType == BodyType2D.Dynamic)
                {
                    total += 0.5 * body.Mass * body.Velocity.MagnitudeSquared;
                    total += 0.5 * body.Inertia * body.AngularVelocity * body.AngularVelocity;
                }
            }
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
}
