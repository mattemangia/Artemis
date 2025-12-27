namespace ArtemisEngine;

public class PhysicsWorld
{
    public Vector2 Gravity { get; set; }
    public List<RigidBody> Bodies { get; private set; }
    public List<Joint> Joints { get; private set; }
    public List<AreaEffector> AreaEffectors { get; private set; }

    // Spatial partitioning for optimization
    private SpatialHashGrid? _spatialGrid;
    public bool UseSpatialPartitioning { get; set; } = true;

    // Collision tracking for events
    private Dictionary<(RigidBody, RigidBody), bool> _previousCollisions;
    private HashSet<(RigidBody, RigidBody)> _currentCollisions;

    // Advanced solver for better stability
    private SequentialImpulseSolver _solver;
    public bool UseAdvancedSolver { get; set; } = true;

    // Events
    public event CollisionEventHandler? OnCollisionEnter;
    public event CollisionEventHandler? OnCollisionStay;
    public event CollisionEventHandler? OnCollisionExit;
    public event CollisionEventHandler? OnTriggerEnter;
    public event CollisionEventHandler? OnTriggerStay;
    public event CollisionEventHandler? OnTriggerExit;

    private const int VelocityIterations = 6;
    private const int PositionIterations = 2;

    public PhysicsWorld(Vector2 gravity)
    {
        Gravity = gravity;
        Bodies = new List<RigidBody>();
        Joints = new List<Joint>();
        AreaEffectors = new List<AreaEffector>();
        _spatialGrid = new SpatialHashGrid(10f);
        _previousCollisions = new Dictionary<(RigidBody, RigidBody), bool>();
        _currentCollisions = new HashSet<(RigidBody, RigidBody)>();
        _solver = new SequentialImpulseSolver();
    }

    public void AddBody(RigidBody body)
    {
        Bodies.Add(body);
    }

    public void RemoveBody(RigidBody body)
    {
        Bodies.Remove(body);

        // Clean up collision tracking
        var keysToRemove = _previousCollisions.Keys
            .Where(k => k.Item1 == body || k.Item2 == body)
            .ToList();

        foreach (var key in keysToRemove)
        {
            _previousCollisions.Remove(key);
        }
    }

    public void AddJoint(Joint joint)
    {
        Joints.Add(joint);
    }

    public void RemoveJoint(Joint joint)
    {
        Joints.Remove(joint);
    }

    public void AddAreaEffector(AreaEffector effector)
    {
        AreaEffectors.Add(effector);
    }

    public void RemoveAreaEffector(AreaEffector effector)
    {
        AreaEffectors.Remove(effector);
    }

    public RaycastHit Raycast(Ray ray, int layerMask = ~0)
    {
        return Raycaster.Raycast(ray, Bodies, layerMask);
    }

    public RaycastHit[] RaycastAll(Ray ray, int layerMask = ~0)
    {
        return Raycaster.RaycastAll(ray, Bodies, layerMask);
    }

    public void Step(float deltaTime)
    {
        _currentCollisions.Clear();

        // Begin solver frame (contact persistence)
        if (UseAdvancedSolver)
        {
            _solver.BeginFrame();
        }

        // Update sleeping states
        foreach (var body in Bodies)
        {
            if (!body.IsStatic && body.CanSleep)
            {
                SleepingSystem.UpdateSleepState(body, deltaTime);
            }
        }

        // Apply area effectors (wind, buoyancy, gravity wells, etc.)
        foreach (var effector in AreaEffectors)
        {
            if (effector.Enabled)
            {
                foreach (var body in Bodies)
                {
                    if (!body.IsStatic && !body.IsSleeping && effector.AffectsBody(body))
                    {
                        effector.ApplyForce(body, deltaTime);
                    }
                }
            }
        }

        // Apply gravity
        foreach (var body in Bodies)
        {
            if (!body.IsStatic && !body.IsKinematic && !body.IsSleeping)
            {
                body.ApplyForce(Gravity * body.Mass * deltaTime);
            }
        }

        // Update velocities and positions
        foreach (var body in Bodies)
        {
            body.Update(deltaTime);
        }

        // Solve joints
        foreach (var joint in Joints)
        {
            if (joint.Enabled)
            {
                joint.Solve(deltaTime);
            }
        }

        // Collision detection and resolution
        if (UseSpatialPartitioning && _spatialGrid != null)
        {
            PerformSpatialCollisionDetection(deltaTime);
        }
        else
        {
            PerformBruteForceCollisionDetection(deltaTime);
        }

        // Process collision events
        ProcessCollisionEvents();
    }

    private void PerformBruteForceCollisionDetection(float deltaTime)
    {
        for (int i = 0; i < Bodies.Count; i++)
        {
            for (int j = i + 1; j < Bodies.Count; j++)
            {
                var bodyA = Bodies[i];
                var bodyB = Bodies[j];

                ProcessCollisionPair(bodyA, bodyB, deltaTime);
            }
        }
    }

    private void PerformSpatialCollisionDetection(float deltaTime)
    {
        // Rebuild spatial grid
        _spatialGrid!.Clear();
        foreach (var body in Bodies)
        {
            _spatialGrid.Insert(body);
        }

        // Check collisions using spatial grid
        var checkedPairs = new HashSet<(RigidBody, RigidBody)>();

        foreach (var body in Bodies)
        {
            if (body.IsSleeping)
                continue;

            var nearby = _spatialGrid.QueryNearby(body);

            foreach (var other in nearby)
            {
                var pair = body.GetHashCode() < other.GetHashCode() ? (body, other) : (other, body);

                if (!checkedPairs.Add(pair))
                    continue;

                ProcessCollisionPair(body, other, deltaTime);
            }
        }
    }

    private void ProcessCollisionPair(RigidBody bodyA, RigidBody bodyB, float deltaTime)
    {
        // Skip if both sleeping
        if (bodyA.IsSleeping && bodyB.IsSleeping)
            return;

        // Check if can collide
        if (!bodyA.CanCollideWith(bodyB))
            return;

        Collision collision;
        bool hasCollision = false;

        // Check for CCD (Continuous Collision Detection) for fast-moving objects
        if ((bodyA.UseCCD || bodyB.UseCCD) && !bodyA.IsStatic && !bodyB.IsStatic)
        {
            // Try swept collision detection
            RigidBody moving = bodyA.UseCCD ? bodyA : bodyB;
            RigidBody target = bodyA.UseCCD ? bodyB : bodyA;

            if (ContinuousCollision.SweptCollision(moving, target, deltaTime, out float toi, out collision))
            {
                // Move body to time of impact
                moving.Position = moving.Position + moving.Velocity * toi * deltaTime;
                hasCollision = true;
            }
        }

        // Fall back to discrete collision detection if CCD didn't find a collision
        if (!hasCollision)
        {
            hasCollision = DetectCollision(bodyA, bodyB, out collision);
        }

        if (hasCollision)
        {
            // Check one-way platform
            if (OneWayPlatformExtensions.ShouldIgnoreCollision(bodyA, bodyB, collision.Normal))
            {
                return; // Pass through one-way platform
            }

            var pair = (bodyA, bodyB);
            _currentCollisions.Add(pair);

            bool isTrigger = bodyA.IsTrigger || bodyB.IsTrigger;

            if (!isTrigger)
            {
                // Physical collision - use advanced solver if enabled
                if (UseAdvancedSolver)
                {
                    _solver.SolveCollision(bodyA, bodyB, collision);
                }
                else
                {
                    ResolveCollision(bodyA, bodyB, collision);
                }

                // Wake up bodies
                bodyA.WakeUp();
                bodyB.WakeUp();

                // Notify collision listeners
                bodyA.CollisionListener?.OnCollisionEnter(bodyB, collision);
                bodyB.CollisionListener?.OnCollisionEnter(bodyA, collision);
            }
            else
            {
                // Trigger collision - no physical resolution
                bodyA.CollisionListener?.OnCollisionEnter(bodyB, collision);
                bodyB.CollisionListener?.OnCollisionEnter(bodyA, collision);
            }
        }
    }

    private void ProcessCollisionEvents()
    {
        // Detect new collisions (Enter)
        foreach (var pair in _currentCollisions)
        {
            if (!_previousCollisions.ContainsKey(pair))
            {
                bool isTrigger = pair.Item1.IsTrigger || pair.Item2.IsTrigger;

                if (isTrigger)
                {
                    OnTriggerEnter?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                }
                else
                {
                    OnCollisionEnter?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                }
            }
            else
            {
                // Collision ongoing (Stay)
                bool isTrigger = pair.Item1.IsTrigger || pair.Item2.IsTrigger;

                if (isTrigger)
                {
                    OnTriggerStay?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                }
                else
                {
                    OnCollisionStay?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                }
            }
        }

        // Detect ended collisions (Exit)
        foreach (var pair in _previousCollisions.Keys)
        {
            if (!_currentCollisions.Contains(pair))
            {
                bool isTrigger = pair.Item1.IsTrigger || pair.Item2.IsTrigger;

                if (isTrigger)
                {
                    OnTriggerExit?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                }
                else
                {
                    OnCollisionExit?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                }

                pair.Item1.CollisionListener?.OnCollisionExit(pair.Item2);
                pair.Item2.CollisionListener?.OnCollisionExit(pair.Item1);
            }
        }

        // Update previous collisions
        _previousCollisions.Clear();
        foreach (var pair in _currentCollisions)
        {
            _previousCollisions[pair] = true;
        }
    }

    private bool DetectCollision(RigidBody bodyA, RigidBody bodyB, out Collision collision)
    {
        collision = new Collision();

        if (bodyA.Shape is CircleShape circleA && bodyB.Shape is CircleShape circleB)
        {
            return CircleVsCircle(bodyA, bodyB, circleA, circleB, out collision);
        }
        else if (bodyA.Shape is BoxShape boxA && bodyB.Shape is BoxShape boxB)
        {
            return BoxVsBox(bodyA, bodyB, boxA, boxB, out collision);
        }
        else if (bodyA.Shape is CircleShape circle && bodyB.Shape is BoxShape box)
        {
            return CircleVsBox(bodyA, bodyB, circle, box, out collision);
        }
        else if (bodyA.Shape is BoxShape box2 && bodyB.Shape is CircleShape circle2)
        {
            bool result = CircleVsBox(bodyB, bodyA, circle2, box2, out collision);
            if (result)
            {
                collision.Normal = -collision.Normal;
            }
            return result;
        }

        return false;
    }

    private bool CircleVsCircle(RigidBody bodyA, RigidBody bodyB, CircleShape circleA, CircleShape circleB, out Collision collision)
    {
        collision = new Collision();

        Vector2 delta = bodyB.Position - bodyA.Position;
        float distanceSquared = delta.LengthSquared;
        float radiusSum = circleA.Radius + circleB.Radius;

        if (distanceSquared >= radiusSum * radiusSum)
            return false;

        float distance = MathF.Sqrt(distanceSquared);

        collision.Normal = distance > 0 ? delta / distance : new Vector2(1, 0);
        collision.Penetration = radiusSum - distance;
        collision.ContactPoint = bodyA.Position + collision.Normal * circleA.Radius;

        return true;
    }

    private bool BoxVsBox(RigidBody bodyA, RigidBody bodyB, BoxShape boxA, BoxShape boxB, out Collision collision)
    {
        // Use SAT (Separating Axis Theorem) for accurate rotated box collision
        return CollisionDetection.BoxVsBoxSAT(bodyA, bodyB, boxA, boxB, out collision);
    }

    private bool CircleVsBox(RigidBody circleBody, RigidBody boxBody, CircleShape circle, BoxShape box, out Collision collision)
    {
        // Use improved circle vs box with rotation support
        return CollisionDetection.CircleVsBoxImproved(circleBody, boxBody, circle, box, out collision);
    }

    private void ResolveCollision(RigidBody bodyA, RigidBody bodyB, Collision collision)
    {
        // Position correction
        float percent = 0.8f;
        float slop = 0.01f;
        Vector2 correction = collision.Normal *
            (Math.Max(collision.Penetration - slop, 0.0f) / (bodyA.InverseMass + bodyB.InverseMass)) * percent;

        if (!bodyA.IsStatic)
            bodyA.Position -= correction * bodyA.InverseMass;
        if (!bodyB.IsStatic)
            bodyB.Position += correction * bodyB.InverseMass;

        // Calculate relative velocity at contact point (includes rotation!)
        Vector2 rA = collision.ContactPoint - bodyA.Position;
        Vector2 rB = collision.ContactPoint - bodyB.Position;

        Vector2 velocityA = bodyA.Velocity + new Vector2(-rA.Y * bodyA.AngularVelocity, rA.X * bodyA.AngularVelocity);
        Vector2 velocityB = bodyB.Velocity + new Vector2(-rB.Y * bodyB.AngularVelocity, rB.X * bodyB.AngularVelocity);
        Vector2 relativeVelocity = velocityB - velocityA;

        float velocityAlongNormal = Vector2.Dot(relativeVelocity, collision.Normal);

        if (velocityAlongNormal > 0)
            return;

        // Calculate impulse with rotation
        float restitution = Math.Min(bodyA.Restitution, bodyB.Restitution);

        float rACrossN = Vector2.Cross(rA, collision.Normal);
        float rBCrossN = Vector2.Cross(rB, collision.Normal);

        float impulseMagnitude = -(1 + restitution) * velocityAlongNormal;
        impulseMagnitude /= bodyA.InverseMass + bodyB.InverseMass +
                           rACrossN * rACrossN * bodyA.InverseInertia +
                           rBCrossN * rBCrossN * bodyB.InverseInertia;

        Vector2 impulse = collision.Normal * impulseMagnitude;

        // Apply impulse at contact point (includes torque!)
        if (!bodyA.IsStatic)
        {
            bodyA.Velocity -= impulse * bodyA.InverseMass;
            bodyA.AngularVelocity -= Vector2.Cross(rA, impulse) * bodyA.InverseInertia;
        }

        if (!bodyB.IsStatic)
        {
            bodyB.Velocity += impulse * bodyB.InverseMass;
            bodyB.AngularVelocity += Vector2.Cross(rB, impulse) * bodyB.InverseInertia;
        }

        // Friction with rotation
        relativeVelocity = velocityB - velocityA;
        Vector2 tangent = relativeVelocity - collision.Normal * Vector2.Dot(relativeVelocity, collision.Normal);

        if (tangent.LengthSquared > 0.0001f)
        {
            tangent = tangent.Normalized;

            float rACrossT = Vector2.Cross(rA, tangent);
            float rBCrossT = Vector2.Cross(rB, tangent);

            float frictionMagnitude = -Vector2.Dot(relativeVelocity, tangent);
            frictionMagnitude /= bodyA.InverseMass + bodyB.InverseMass +
                                rACrossT * rACrossT * bodyA.InverseInertia +
                                rBCrossT * rBCrossT * bodyB.InverseInertia;

            float mu = MathF.Sqrt(bodyA.Friction * bodyA.Friction + bodyB.Friction * bodyB.Friction);
            Vector2 frictionImpulse = tangent * Math.Clamp(frictionMagnitude, -impulseMagnitude * mu, impulseMagnitude * mu);

            if (!bodyA.IsStatic)
            {
                bodyA.Velocity -= frictionImpulse * bodyA.InverseMass;
                bodyA.AngularVelocity -= Vector2.Cross(rA, frictionImpulse) * bodyA.InverseInertia;
            }

            if (!bodyB.IsStatic)
            {
                bodyB.Velocity += frictionImpulse * bodyB.InverseMass;
                bodyB.AngularVelocity += Vector2.Cross(rB, frictionImpulse) * bodyB.InverseInertia;
            }
        }
    }
}

public struct Collision
{
    public Vector2 Normal;
    public float Penetration;
    public float PenetrationDepth => Penetration; // Alias for compatibility
    public Vector2 ContactPoint;
}
