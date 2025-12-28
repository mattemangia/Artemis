using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace ArtemisEngine;

/// <summary>
/// High-performance 2D physics world with multi-threading support.
/// </summary>
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
    private ConcurrentDictionary<(RigidBody, RigidBody), bool> _currentCollisions;

    // Advanced solver for better stability
    private SequentialImpulseSolver _solver;
    public bool UseAdvancedSolver { get; set; } = true;

    // Thread-safe collision accumulator
    private ConcurrentBag<(RigidBody, RigidBody, Collision)> _detectedCollisions;

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
        _currentCollisions = new ConcurrentDictionary<(RigidBody, RigidBody), bool>();
        _detectedCollisions = new ConcurrentBag<(RigidBody, RigidBody, Collision)>();
        _solver = new SequentialImpulseSolver();
    }

    public void AddBody(RigidBody body)
    {
        Bodies.Add(body);
    }

    public void RemoveBody(RigidBody body)
    {
        Bodies.Remove(body);

        var keysToRemove = _previousCollisions.Keys
            .Where(k => k.Item1 == body || k.Item2 == body)
            .ToList();

        foreach (var key in keysToRemove)
        {
            _previousCollisions.Remove(key);
        }
    }

    public void AddJoint(Joint joint) => Joints.Add(joint);
    public void RemoveJoint(Joint joint) => Joints.Remove(joint);
    public void AddAreaEffector(AreaEffector effector) => AreaEffectors.Add(effector);
    public void RemoveAreaEffector(AreaEffector effector) => AreaEffectors.Remove(effector);

    public RaycastHit Raycast(Ray ray, int layerMask = ~0)
        => Raycaster.Raycast(ray, Bodies, layerMask);

    public RaycastHit[] RaycastAll(Ray ray, int layerMask = ~0)
        => Raycaster.RaycastAll(ray, Bodies, layerMask);

    public void Step(float deltaTime)
    {
        _currentCollisions.Clear();
        while (_detectedCollisions.TryTake(out _)) { }

        if (UseAdvancedSolver)
        {
            _solver.BeginFrame();
        }

        // Parallel: Update sleeping states
        var activeBodies = Bodies.Where(b => !b.IsStatic && b.CanSleep).ToArray();
        Parallel.ForEach(activeBodies, body =>
        {
            SleepingSystem.UpdateSleepState(body, deltaTime);
        });

        // Parallel: Apply area effectors
        var awakeBodies = Bodies.Where(b => !b.IsStatic && !b.IsSleeping).ToArray();
        var enabledEffectors = AreaEffectors.Where(e => e.Enabled).ToArray();

        Parallel.ForEach(awakeBodies, body =>
        {
            foreach (var effector in enabledEffectors)
            {
                if (effector.AffectsBody(body))
                {
                    effector.ApplyForce(body, deltaTime);
                }
            }
        });

        // Parallel: Apply gravity and integrate
        var dynamicBodies = Bodies.Where(b => !b.IsStatic && !b.IsKinematic && !b.IsSleeping).ToArray();
        var gravity = Gravity;

        Parallel.ForEach(dynamicBodies, body =>
        {
            body.ApplyForce(gravity * body.Mass * deltaTime);
        });

        // Parallel: Update velocities and positions
        Parallel.ForEach(Bodies, body =>
        {
            body.Update(deltaTime);
        });

        // Sequential: Solve joints (order matters for constraints)
        foreach (var joint in Joints)
        {
            if (joint.Enabled)
            {
                joint.Solve(deltaTime);
            }
        }

        // Collision detection (parallel detection, sequential resolution)
        if (UseSpatialPartitioning && _spatialGrid != null)
        {
            PerformSpatialCollisionDetectionParallel(deltaTime);
        }
        else
        {
            PerformBruteForceCollisionDetectionParallel(deltaTime);
        }

        // Process collision events
        ProcessCollisionEvents();
    }

    private void PerformBruteForceCollisionDetectionParallel(float deltaTime)
    {
        var bodiesArray = Bodies.ToArray();
        int count = bodiesArray.Length;

        // Parallel collision detection
        Parallel.For(0, count, i =>
        {
            for (int j = i + 1; j < count; j++)
            {
                DetectCollisionPair(bodiesArray[i], bodiesArray[j], deltaTime);
            }
        });

        // Sequential collision resolution (for determinism)
        ResolveDetectedCollisions(deltaTime);
    }

    private void PerformSpatialCollisionDetectionParallel(float deltaTime)
    {
        // Rebuild spatial grid (sequential - grid not thread-safe for writes)
        _spatialGrid!.Clear();
        foreach (var body in Bodies)
        {
            _spatialGrid.Insert(body);
        }

        // Get awake bodies for parallel processing
        var awakeBodies = Bodies.Where(b => !b.IsSleeping).ToArray();
        var checkedPairs = new ConcurrentDictionary<(RigidBody, RigidBody), bool>();

        // Parallel collision detection
        Parallel.ForEach(awakeBodies, body =>
        {
            var nearby = _spatialGrid.QueryNearby(body);

            foreach (var other in nearby)
            {
                var pair = body.GetHashCode() < other.GetHashCode() ? (body, other) : (other, body);

                if (!checkedPairs.TryAdd(pair, true))
                    continue;

                DetectCollisionPair(body, other, deltaTime);
            }
        });

        // Sequential collision resolution
        ResolveDetectedCollisions(deltaTime);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DetectCollisionPair(RigidBody bodyA, RigidBody bodyB, float deltaTime)
    {
        if (bodyA.IsSleeping && bodyB.IsSleeping) return;
        if (!bodyA.CanCollideWith(bodyB)) return;

        Collision collision = default;
        bool hasCollision = false;

        // CCD for fast-moving objects
        if ((bodyA.UseCCD || bodyB.UseCCD) && !bodyA.IsStatic && !bodyB.IsStatic)
        {
            RigidBody moving = bodyA.UseCCD ? bodyA : bodyB;
            RigidBody target = bodyA.UseCCD ? bodyB : bodyA;

            if (ContinuousCollision.SweptCollision(moving, target, deltaTime, out float toi, out collision))
            {
                hasCollision = true;
            }
        }

        if (!hasCollision)
        {
            hasCollision = DetectCollision(bodyA, bodyB, out collision);
        }

        if (hasCollision)
        {
            if (OneWayPlatformExtensions.ShouldIgnoreCollision(bodyA, bodyB, collision.Normal))
                return;

            _detectedCollisions.Add((bodyA, bodyB, collision));
        }
    }

    private void ResolveDetectedCollisions(float deltaTime)
    {
        foreach (var (bodyA, bodyB, collision) in _detectedCollisions)
        {
            var pair = (bodyA, bodyB);
            _currentCollisions.TryAdd(pair, true);

            bool isTrigger = bodyA.IsTrigger || bodyB.IsTrigger;

            if (!isTrigger)
            {
                // Handle CCD position update
                if ((bodyA.UseCCD || bodyB.UseCCD) && !bodyA.IsStatic && !bodyB.IsStatic)
                {
                    RigidBody moving = bodyA.UseCCD ? bodyA : bodyB;
                    // Position already interpolated during detection
                }

                if (UseAdvancedSolver)
                {
                    _solver.SolveCollision(bodyA, bodyB, collision);
                }
                else
                {
                    ResolveCollision(bodyA, bodyB, collision);
                }

                bodyA.WakeUp();
                bodyB.WakeUp();

                bodyA.CollisionListener?.OnCollisionEnter(bodyB, collision);
                bodyB.CollisionListener?.OnCollisionEnter(bodyA, collision);
            }
            else
            {
                bodyA.CollisionListener?.OnCollisionEnter(bodyB, collision);
                bodyB.CollisionListener?.OnCollisionEnter(bodyA, collision);
            }
        }
    }

    private void ProcessCollisionEvents()
    {
        var currentKeys = _currentCollisions.Keys.ToHashSet();

        // Detect new collisions (Enter) and ongoing (Stay)
        foreach (var pair in currentKeys)
        {
            if (!_previousCollisions.ContainsKey(pair))
            {
                bool isTrigger = pair.Item1.IsTrigger || pair.Item2.IsTrigger;

                if (isTrigger)
                    OnTriggerEnter?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                else
                    OnCollisionEnter?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
            }
            else
            {
                bool isTrigger = pair.Item1.IsTrigger || pair.Item2.IsTrigger;

                if (isTrigger)
                    OnTriggerStay?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                else
                    OnCollisionStay?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
            }
        }

        // Detect ended collisions (Exit)
        foreach (var pair in _previousCollisions.Keys)
        {
            if (!currentKeys.Contains(pair))
            {
                bool isTrigger = pair.Item1.IsTrigger || pair.Item2.IsTrigger;

                if (isTrigger)
                    OnTriggerExit?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));
                else
                    OnCollisionExit?.Invoke(this, new CollisionEventArgs(pair.Item1, pair.Item2, new Collision()));

                pair.Item1.CollisionListener?.OnCollisionExit(pair.Item2);
                pair.Item2.CollisionListener?.OnCollisionExit(pair.Item1);
            }
        }

        // Update previous collisions
        _previousCollisions.Clear();
        foreach (var pair in currentKeys)
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

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private bool CircleVsCircle(RigidBody bodyA, RigidBody bodyB, CircleShape circleA, CircleShape circleB, out Collision collision)
    {
        collision = new Collision();

        float dx = bodyB.Position.X - bodyA.Position.X;
        float dy = bodyB.Position.Y - bodyA.Position.Y;
        float distanceSquared = dx * dx + dy * dy;
        float radiusSum = circleA.Radius + circleB.Radius;

        if (distanceSquared >= radiusSum * radiusSum)
            return false;

        float distance = MathF.Sqrt(distanceSquared);

        if (distance > 0)
        {
            float invDist = 1f / distance;
            collision.Normal = new Vector2(dx * invDist, dy * invDist);
        }
        else
        {
            collision.Normal = new Vector2(1, 0);
        }

        collision.Penetration = radiusSum - distance;
        collision.ContactPoint = bodyA.Position + collision.Normal * circleA.Radius;

        return true;
    }

    private bool BoxVsBox(RigidBody bodyA, RigidBody bodyB, BoxShape boxA, BoxShape boxB, out Collision collision)
    {
        return CollisionDetection.BoxVsBoxSAT(bodyA, bodyB, boxA, boxB, out collision);
    }

    private bool CircleVsBox(RigidBody circleBody, RigidBody boxBody, CircleShape circle, BoxShape box, out Collision collision)
    {
        return CollisionDetection.CircleVsBoxImproved(circleBody, boxBody, circle, box, out collision);
    }

    private void ResolveCollision(RigidBody bodyA, RigidBody bodyB, Collision collision)
    {
        float percent = 0.8f;
        float slop = 0.01f;
        float correction = Math.Max(collision.Penetration - slop, 0.0f) / (bodyA.InverseMass + bodyB.InverseMass) * percent;
        Vector2 correctionVec = collision.Normal * correction;

        if (!bodyA.IsStatic)
            bodyA.Position -= correctionVec * bodyA.InverseMass;
        if (!bodyB.IsStatic)
            bodyB.Position += correctionVec * bodyB.InverseMass;

        Vector2 rA = collision.ContactPoint - bodyA.Position;
        Vector2 rB = collision.ContactPoint - bodyB.Position;

        Vector2 velocityA = bodyA.Velocity + new Vector2(-rA.Y * bodyA.AngularVelocity, rA.X * bodyA.AngularVelocity);
        Vector2 velocityB = bodyB.Velocity + new Vector2(-rB.Y * bodyB.AngularVelocity, rB.X * bodyB.AngularVelocity);
        Vector2 relativeVelocity = velocityB - velocityA;

        float velocityAlongNormal = Vector2.Dot(relativeVelocity, collision.Normal);

        if (velocityAlongNormal > 0)
            return;

        float restitution = Math.Min(bodyA.Restitution, bodyB.Restitution);

        float rACrossN = Vector2.Cross(rA, collision.Normal);
        float rBCrossN = Vector2.Cross(rB, collision.Normal);

        float impulseMagnitude = -(1 + restitution) * velocityAlongNormal;
        impulseMagnitude /= bodyA.InverseMass + bodyB.InverseMass +
                           rACrossN * rACrossN * bodyA.InverseInertia +
                           rBCrossN * rBCrossN * bodyB.InverseInertia;

        Vector2 impulse = collision.Normal * impulseMagnitude;

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

        // Friction
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
    public float PenetrationDepth => Penetration;
    public Vector2 ContactPoint;
}
