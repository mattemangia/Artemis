namespace ArtemisEngine;

/// <summary>
/// Contact manifold - stores persistent contact information between frames
/// Used for warm starting to improve solver convergence
/// </summary>
public class ContactManifold
{
    public RigidBody BodyA { get; set; }
    public RigidBody BodyB { get; set; }
    public List<ContactPoint> ContactPoints { get; set; } = new();
    public int FrameCount { get; set; } = 0; // How many frames this contact has persisted

    public class ContactPoint
    {
        public Vector2 LocalPointA { get; set; } // Contact point in body A's local space
        public Vector2 LocalPointB { get; set; } // Contact point in body B's local space
        public float NormalImpulse { get; set; } // Accumulated normal impulse
        public float TangentImpulse { get; set; } // Accumulated friction impulse
        public float Separation { get; set; } // Penetration depth (negative if penetrating)
        public Vector2 Normal { get; set; }
    }

    public string GetKey()
    {
        // Create unique key for this contact pair (order-independent)
        int hashA = BodyA.GetHashCode();
        int hashB = BodyB.GetHashCode();
        return hashA < hashB ? $"{hashA}_{hashB}" : $"{hashB}_{hashA}";
    }
}

/// <summary>
/// Contact persistence - maintains contact information across frames
/// Enables warm starting for faster convergence
/// </summary>
public class ContactPersistence
{
    private Dictionary<string, ContactManifold> _manifolds = new();
    private const int MaxFrameAge = 3; // Remove contacts that haven't been updated in 3 frames
    private const float ContactMatchThreshold = 0.1f; // Distance to consider contacts the same

    public void BeginFrame()
    {
        // Age out old contacts
        var keysToRemove = new List<string>();
        foreach (var pair in _manifolds)
        {
            pair.Value.FrameCount++;
            if (pair.Value.FrameCount > MaxFrameAge)
            {
                keysToRemove.Add(pair.Key);
            }
        }

        foreach (var key in keysToRemove)
        {
            _manifolds.Remove(key);
        }
    }

    public ContactManifold GetOrCreateManifold(RigidBody bodyA, RigidBody bodyB)
    {
        var tempManifold = new ContactManifold { BodyA = bodyA, BodyB = bodyB };
        string key = tempManifold.GetKey();

        if (!_manifolds.TryGetValue(key, out var manifold))
        {
            manifold = new ContactManifold { BodyA = bodyA, BodyB = bodyB };
            _manifolds[key] = manifold;
        }

        // Reset frame count - this contact is active
        manifold.FrameCount = 0;
        return manifold;
    }

    /// <summary>
    /// Match new contact points with existing ones for warm starting
    /// </summary>
    public void UpdateManifold(ContactManifold manifold, Vector2 worldContactPoint, Vector2 normal, float separation)
    {
        // Convert to local space
        Vector2 localA = WorldToLocal(manifold.BodyA, worldContactPoint);
        Vector2 localB = WorldToLocal(manifold.BodyB, worldContactPoint);

        // Try to find matching existing contact point
        ContactManifold.ContactPoint? matched = null;
        foreach (var cp in manifold.ContactPoints)
        {
            Vector2 worldA = LocalToWorld(manifold.BodyA, cp.LocalPointA);
            float dist = (worldA - worldContactPoint).Length;

            if (dist < ContactMatchThreshold)
            {
                matched = cp;
                break;
            }
        }

        if (matched != null)
        {
            // Update existing contact (keep accumulated impulses for warm starting!)
            matched.LocalPointA = localA;
            matched.LocalPointB = localB;
            matched.Separation = separation;
            matched.Normal = normal;
        }
        else
        {
            // New contact point
            manifold.ContactPoints.Add(new ContactManifold.ContactPoint
            {
                LocalPointA = localA,
                LocalPointB = localB,
                NormalImpulse = 0,
                TangentImpulse = 0,
                Separation = separation,
                Normal = normal
            });
        }
    }

    private Vector2 WorldToLocal(RigidBody body, Vector2 worldPoint)
    {
        Vector2 relative = worldPoint - body.Position;
        float cos = MathF.Cos(-body.Rotation);
        float sin = MathF.Sin(-body.Rotation);
        return new Vector2(
            relative.X * cos - relative.Y * sin,
            relative.X * sin + relative.Y * cos
        );
    }

    private Vector2 LocalToWorld(RigidBody body, Vector2 localPoint)
    {
        float cos = MathF.Cos(body.Rotation);
        float sin = MathF.Sin(body.Rotation);
        return new Vector2(
            localPoint.X * cos - localPoint.Y * sin + body.Position.X,
            localPoint.X * sin + localPoint.Y * cos + body.Position.Y
        );
    }

    public void Clear()
    {
        _manifolds.Clear();
    }
}

/// <summary>
/// Iterative constraint solver using Sequential Impulse method
/// Much more stable than single-pass impulse resolution
/// </summary>
public class SequentialImpulseSolver
{
    public int VelocityIterations { get; set; } = 8;
    public int PositionIterations { get; set; } = 3;
    public float Baumgarte { get; set; } = 0.2f; // Position correction factor
    public float SlowdownFactor { get; set; } = 1.0f; // For penetration resolution

    private ContactPersistence _persistence = new();

    public void BeginFrame()
    {
        _persistence.BeginFrame();
    }

    /// <summary>
    /// Solve collision using iterative Sequential Impulse method with warm starting
    /// </summary>
    public void SolveCollision(RigidBody bodyA, RigidBody bodyB, Collision collision)
    {
        if (bodyA.IsStatic && bodyB.IsStatic)
            return;

        Vector2 normal = collision.Normal;
        Vector2 contactPoint = collision.ContactPoint;
        float penetration = collision.PenetrationDepth;

        // Get or create contact manifold for persistence
        var manifold = _persistence.GetOrCreateManifold(bodyA, bodyB);
        _persistence.UpdateManifold(manifold, contactPoint, normal, -penetration);

        // Calculate relative position from center of mass
        Vector2 rA = contactPoint - bodyA.Position;
        Vector2 rB = contactPoint - bodyB.Position;

        // Calculate effective mass for normal constraint
        float rnA = Vector2.Cross(rA, normal);
        float rnB = Vector2.Cross(rB, normal);
        float effectiveMassNormal = bodyA.InverseMass + bodyB.InverseMass +
                                   rnA * rnA * bodyA.InverseInertia +
                                   rnB * rnB * bodyB.InverseInertia;

        if (effectiveMassNormal < 0.0001f)
            return;

        effectiveMassNormal = 1.0f / effectiveMassNormal;

        // Calculate tangent (friction direction)
        Vector2 tangent = new Vector2(-normal.Y, normal.X);

        // Calculate effective mass for tangent constraint
        float rtA = Vector2.Cross(rA, tangent);
        float rtB = Vector2.Cross(rB, tangent);
        float effectiveMassTangent = bodyA.InverseMass + bodyB.InverseMass +
                                    rtA * rtA * bodyA.InverseInertia +
                                    rtB * rtB * bodyB.InverseInertia;

        effectiveMassTangent = effectiveMassTangent < 0.0001f ? 0 : 1.0f / effectiveMassTangent;

        // Combined restitution and friction
        float restitution = MathF.Max(bodyA.Restitution, bodyB.Restitution);
        float friction = MathF.Sqrt(bodyA.Friction * bodyB.Friction); // Geometric mean

        // Get cached impulses for warm starting
        float accumulatedNormalImpulse = 0;
        float accumulatedTangentImpulse = 0;

        if (manifold.ContactPoints.Count > 0)
        {
            var cp = manifold.ContactPoints[0];
            accumulatedNormalImpulse = cp.NormalImpulse;
            accumulatedTangentImpulse = cp.TangentImpulse;
        }

        // Warm start - apply cached impulses
        if (accumulatedNormalImpulse > 0 || MathF.Abs(accumulatedTangentImpulse) > 0)
        {
            Vector2 impulse = normal * accumulatedNormalImpulse + tangent * accumulatedTangentImpulse;

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
        }

        // Velocity iterations (resolve velocity constraints)
        for (int i = 0; i < VelocityIterations; i++)
        {
            // Calculate relative velocity at contact point
            Vector2 vA = bodyA.Velocity + new Vector2(-rA.Y * bodyA.AngularVelocity, rA.X * bodyA.AngularVelocity);
            Vector2 vB = bodyB.Velocity + new Vector2(-rB.Y * bodyB.AngularVelocity, rB.X * bodyB.AngularVelocity);
            Vector2 relativeVelocity = vB - vA;

            // Normal impulse (collision response)
            float normalVelocity = Vector2.Dot(relativeVelocity, normal);

            // Apply restitution only if objects are approaching
            float restitutionBias = normalVelocity < -1.0f ? -restitution * normalVelocity : 0;

            float normalImpulseMagnitude = (-normalVelocity + restitutionBias) * effectiveMassNormal;

            // Clamp accumulated impulse (can only push apart, not pull together)
            float newAccumulatedNormal = MathF.Max(accumulatedNormalImpulse + normalImpulseMagnitude, 0);
            normalImpulseMagnitude = newAccumulatedNormal - accumulatedNormalImpulse;
            accumulatedNormalImpulse = newAccumulatedNormal;

            Vector2 normalImpulse = normal * normalImpulseMagnitude;

            // Apply normal impulse
            if (!bodyA.IsStatic)
            {
                bodyA.Velocity -= normalImpulse * bodyA.InverseMass;
                bodyA.AngularVelocity -= Vector2.Cross(rA, normalImpulse) * bodyA.InverseInertia;
            }

            if (!bodyB.IsStatic)
            {
                bodyB.Velocity += normalImpulse * bodyB.InverseMass;
                bodyB.AngularVelocity += Vector2.Cross(rB, normalImpulse) * bodyB.InverseInertia;
            }

            // Friction impulse (Coulomb friction)
            if (effectiveMassTangent > 0 && friction > 0)
            {
                // Recalculate relative velocity after normal impulse
                vA = bodyA.Velocity + new Vector2(-rA.Y * bodyA.AngularVelocity, rA.X * bodyA.AngularVelocity);
                vB = bodyB.Velocity + new Vector2(-rB.Y * bodyB.AngularVelocity, rB.X * bodyB.AngularVelocity);
                relativeVelocity = vB - vA;

                float tangentVelocity = Vector2.Dot(relativeVelocity, tangent);
                float tangentImpulseMagnitude = -tangentVelocity * effectiveMassTangent;

                // Coulomb friction: clamp friction impulse to friction cone
                float maxFriction = friction * accumulatedNormalImpulse;
                float newAccumulatedTangent = MathF.Max(-maxFriction, MathF.Min(accumulatedTangentImpulse + tangentImpulseMagnitude, maxFriction));
                tangentImpulseMagnitude = newAccumulatedTangent - accumulatedTangentImpulse;
                accumulatedTangentImpulse = newAccumulatedTangent;

                Vector2 frictionImpulse = tangent * tangentImpulseMagnitude;

                // Apply friction impulse
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

        // Store accumulated impulses for next frame (warm starting)
        if (manifold.ContactPoints.Count > 0)
        {
            manifold.ContactPoints[0].NormalImpulse = accumulatedNormalImpulse;
            manifold.ContactPoints[0].TangentImpulse = accumulatedTangentImpulse;
        }

        // Position iterations (resolve penetration)
        for (int i = 0; i < PositionIterations; i++)
        {
            if (penetration <= 0)
                break;

            // Calculate position correction
            float correction = MathF.Max(penetration - 0.01f, 0) * Baumgarte;
            Vector2 positionImpulse = normal * correction;

            // Apply position correction
            if (!bodyA.IsStatic)
            {
                float massRatioA = bodyA.InverseMass / (bodyA.InverseMass + bodyB.InverseMass);
                bodyA.Position -= positionImpulse * massRatioA * SlowdownFactor;
            }

            if (!bodyB.IsStatic)
            {
                float massRatioB = bodyB.InverseMass / (bodyA.InverseMass + bodyB.InverseMass);
                bodyB.Position += positionImpulse * massRatioB * SlowdownFactor;
            }

            penetration *= 0.5f; // Reduce penetration for next iteration
        }
    }

    public void Clear()
    {
        _persistence.Clear();
    }
}
