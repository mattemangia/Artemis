namespace ArtemisEngine;

/// <summary>
/// Base class for physics joints/constraints that connect two bodies
/// </summary>
public abstract class Joint
{
    public RigidBody BodyA { get; set; }
    public RigidBody BodyB { get; set; }
    public bool Enabled { get; set; } = true;
    public float Stiffness { get; set; } = 1.0f;

    protected Joint(RigidBody bodyA, RigidBody bodyB)
    {
        BodyA = bodyA;
        BodyB = bodyB;
    }

    public abstract void Solve(float deltaTime);
}

/// <summary>
/// Distance joint - maintains a fixed distance between two bodies
/// </summary>
public class DistanceJoint : Joint
{
    public float Distance { get; set; }
    public Vector2 LocalAnchorA { get; set; }
    public Vector2 LocalAnchorB { get; set; }

    public DistanceJoint(RigidBody bodyA, RigidBody bodyB, float distance)
        : base(bodyA, bodyB)
    {
        Distance = distance;
        LocalAnchorA = Vector2.Zero;
        LocalAnchorB = Vector2.Zero;
    }

    public DistanceJoint(RigidBody bodyA, RigidBody bodyB, Vector2 localAnchorA, Vector2 localAnchorB)
        : base(bodyA, bodyB)
    {
        LocalAnchorA = localAnchorA;
        LocalAnchorB = localAnchorB;
        Vector2 worldAnchorA = bodyA.Position + localAnchorA;
        Vector2 worldAnchorB = bodyB.Position + localAnchorB;
        Distance = (worldAnchorB - worldAnchorA).Length;
    }

    public override void Solve(float deltaTime)
    {
        if (!Enabled) return;

        Vector2 anchorA = BodyA.Position + LocalAnchorA;
        Vector2 anchorB = BodyB.Position + LocalAnchorB;

        Vector2 delta = anchorB - anchorA;
        float currentDistance = delta.Length;

        if (currentDistance < 0.0001f)
            return;

        float error = currentDistance - Distance;
        Vector2 direction = delta / currentDistance;

        float totalInverseMass = BodyA.InverseMass + BodyB.InverseMass;
        if (totalInverseMass < 0.0001f)
            return;

        Vector2 correction = direction * (error / totalInverseMass) * Stiffness;

        if (!BodyA.IsStatic)
            BodyA.Position += correction * BodyA.InverseMass;
        if (!BodyB.IsStatic)
            BodyB.Position -= correction * BodyB.InverseMass;
    }
}

/// <summary>
/// Spring joint - acts like a spring between two bodies
/// </summary>
public class SpringJoint : Joint
{
    public float RestLength { get; set; }
    public float SpringConstant { get; set; }
    public float Damping { get; set; }
    public Vector2 LocalAnchorA { get; set; }
    public Vector2 LocalAnchorB { get; set; }

    public SpringJoint(RigidBody bodyA, RigidBody bodyB, float restLength, float springConstant, float damping = 0.1f)
        : base(bodyA, bodyB)
    {
        RestLength = restLength;
        SpringConstant = springConstant;
        Damping = damping;
        LocalAnchorA = Vector2.Zero;
        LocalAnchorB = Vector2.Zero;
    }

    public override void Solve(float deltaTime)
    {
        if (!Enabled) return;

        Vector2 anchorA = BodyA.Position + LocalAnchorA;
        Vector2 anchorB = BodyB.Position + LocalAnchorB;

        Vector2 delta = anchorB - anchorA;
        float distance = delta.Length;

        if (distance < 0.0001f)
            return;

        Vector2 direction = delta / distance;

        // Spring force: F = -k * x
        float displacement = distance - RestLength;
        Vector2 springForce = direction * (SpringConstant * displacement);

        // Damping force
        Vector2 relativeVelocity = BodyB.Velocity - BodyA.Velocity;
        Vector2 dampingForce = relativeVelocity * Damping;

        Vector2 totalForce = springForce + dampingForce;

        if (!BodyA.IsStatic)
            BodyA.ApplyForce(totalForce * deltaTime);
        if (!BodyB.IsStatic)
            BodyB.ApplyForce(-totalForce * deltaTime);
    }
}

/// <summary>
/// Revolute joint - allows rotation around a fixed point
/// </summary>
public class RevoluteJoint : Joint
{
    public Vector2 Anchor { get; set; }
    public Vector2 LocalAnchorA { get; set; }
    public Vector2 LocalAnchorB { get; set; }

    public RevoluteJoint(RigidBody bodyA, RigidBody bodyB, Vector2 worldAnchor)
        : base(bodyA, bodyB)
    {
        Anchor = worldAnchor;
        LocalAnchorA = worldAnchor - bodyA.Position;
        LocalAnchorB = worldAnchor - bodyB.Position;
    }

    public override void Solve(float deltaTime)
    {
        if (!Enabled) return;

        Vector2 anchorA = BodyA.Position + LocalAnchorA;
        Vector2 anchorB = BodyB.Position + LocalAnchorB;

        Vector2 delta = anchorB - anchorA;

        float totalInverseMass = BodyA.InverseMass + BodyB.InverseMass;
        if (totalInverseMass < 0.0001f)
            return;

        Vector2 correction = delta / totalInverseMass * Stiffness;

        if (!BodyA.IsStatic)
            BodyA.Position += correction * BodyA.InverseMass;
        if (!BodyB.IsStatic)
            BodyB.Position -= correction * BodyB.InverseMass;
    }
}

/// <summary>
/// Rope joint - constrains maximum distance between bodies
/// </summary>
public class RopeJoint : Joint
{
    public float MaxLength { get; set; }
    public Vector2 LocalAnchorA { get; set; }
    public Vector2 LocalAnchorB { get; set; }

    public RopeJoint(RigidBody bodyA, RigidBody bodyB, float maxLength)
        : base(bodyA, bodyB)
    {
        MaxLength = maxLength;
        LocalAnchorA = Vector2.Zero;
        LocalAnchorB = Vector2.Zero;
    }

    public override void Solve(float deltaTime)
    {
        if (!Enabled) return;

        Vector2 anchorA = BodyA.Position + LocalAnchorA;
        Vector2 anchorB = BodyB.Position + LocalAnchorB;

        Vector2 delta = anchorB - anchorA;
        float currentLength = delta.Length;

        if (currentLength <= MaxLength)
            return;

        float excess = currentLength - MaxLength;
        Vector2 direction = delta / currentLength;

        float totalInverseMass = BodyA.InverseMass + BodyB.InverseMass;
        if (totalInverseMass < 0.0001f)
            return;

        Vector2 correction = direction * (excess / totalInverseMass) * Stiffness;

        if (!BodyA.IsStatic)
            BodyA.Position += correction * BodyA.InverseMass;
        if (!BodyB.IsStatic)
            BodyB.Position -= correction * BodyB.InverseMass;
    }
}
