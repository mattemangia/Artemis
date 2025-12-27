namespace ArtemisEngine;

/// <summary>
/// Base class for area effectors that apply forces to bodies in a region
/// </summary>
public abstract class AreaEffector
{
    public Vector2 Position { get; set; }
    public float Radius { get; set; }
    public bool Enabled { get; set; } = true;
    public int AffectedLayers { get; set; } = CollisionLayers.Everything;

    protected AreaEffector(Vector2 position, float radius)
    {
        Position = position;
        Radius = radius;
    }

    public abstract void ApplyEffect(RigidBody body, float deltaTime);

    protected bool IsBodyAffected(RigidBody body)
    {
        if (!Enabled || body.IsStatic)
            return false;

        // Check layer
        if ((body.CollisionLayer & AffectedLayers) == 0)
            return false;

        // Check distance
        float distanceSquared = (body.Position - Position).LengthSquared;
        return distanceSquared <= Radius * Radius;
    }

    protected float GetFalloff(RigidBody body)
    {
        float distance = (body.Position - Position).Length;
        return 1.0f - (distance / Radius);
    }
}

/// <summary>
/// Applies constant force in a direction (wind, current)
/// </summary>
public class DirectionalForceEffector : AreaEffector
{
    public Vector2 Force { get; set; }
    public bool UseFalloff { get; set; } = true;

    public DirectionalForceEffector(Vector2 position, float radius, Vector2 force)
        : base(position, radius)
    {
        Force = force;
    }

    public override void ApplyEffect(RigidBody body, float deltaTime)
    {
        if (!IsBodyAffected(body))
            return;

        Vector2 appliedForce = Force;

        if (UseFalloff)
        {
            appliedForce *= GetFalloff(body);
        }

        body.ApplyForce(appliedForce * deltaTime);
    }
}

/// <summary>
/// Applies force toward or away from center (gravity well, explosion)
/// </summary>
public class RadialForceEffector : AreaEffector
{
    public float Strength { get; set; }
    public bool IsAttraction { get; set; } = true; // false for repulsion

    public RadialForceEffector(Vector2 position, float radius, float strength, bool attraction = true)
        : base(position, radius)
    {
        Strength = strength;
        IsAttraction = attraction;
    }

    public override void ApplyEffect(RigidBody body, float deltaTime)
    {
        if (!IsBodyAffected(body))
            return;

        Vector2 direction = Position - body.Position;
        float distance = direction.Length;

        if (distance < 0.001f)
            return;

        direction = direction / distance; // Normalize

        if (!IsAttraction)
            direction = -direction;

        // Inverse square falloff like real gravity
        float forceMagnitude = Strength / (distance * distance + 1.0f);
        Vector2 force = direction * forceMagnitude;

        body.ApplyForce(force * deltaTime);
    }
}

/// <summary>
/// Applies buoyancy force (water, lava)
/// </summary>
public class BuoyancyEffector : AreaEffector
{
    public float Density { get; set; } = 1.0f; // Water = 1.0
    public float LinearDrag { get; set; } = 0.5f;
    public float AngularDrag { get; set; } = 0.5f;
    public Vector2 FlowVelocity { get; set; } = Vector2.Zero;
    public float SurfaceY { get; set; }

    public BuoyancyEffector(Vector2 position, float radius, float surfaceY)
        : base(position, radius)
    {
        SurfaceY = surfaceY;
    }

    public override void ApplyEffect(RigidBody body, float deltaTime)
    {
        if (!IsBodyAffected(body))
            return;

        // Simple buoyancy: upward force proportional to submerged volume
        float submergedDepth = SurfaceY - body.Position.Y;

        if (submergedDepth > 0)
        {
            // Buoyant force (simplified)
            float buoyancy = submergedDepth * Density * 10f; // 10 = gravity constant
            body.ApplyForce(new Vector2(0, buoyancy) * deltaTime);

            // Drag forces
            Vector2 dragForce = -body.Velocity * LinearDrag * body.Mass;
            body.ApplyForce(dragForce * deltaTime);

            float angularDrag = -body.AngularVelocity * AngularDrag * body.Inertia;
            body.ApplyTorque(angularDrag * deltaTime);

            // Flow/current
            if (FlowVelocity.LengthSquared > 0)
            {
                body.ApplyForce(FlowVelocity * Density * deltaTime);
            }
        }
    }
}

/// <summary>
/// Applies torque to make objects spin (tornado, vortex)
/// </summary>
public class VortexEffector : AreaEffector
{
    public float AngularForce { get; set; }
    public bool Clockwise { get; set; } = true;

    public VortexEffector(Vector2 position, float radius, float angularForce)
        : base(position, radius)
    {
        AngularForce = angularForce;
    }

    public override void ApplyEffect(RigidBody body, float deltaTime)
    {
        if (!IsBodyAffected(body))
            return;

        Vector2 toCenter = Position - body.Position;
        float distance = toCenter.Length;

        if (distance < 0.001f)
            return;

        // Tangential force (perpendicular to radius)
        Vector2 tangent = new Vector2(-toCenter.Y, toCenter.X).Normalized;

        if (!Clockwise)
            tangent = -tangent;

        float forceMagnitude = AngularForce * GetFalloff(body);
        body.ApplyForce(tangent * forceMagnitude * deltaTime);

        // Also apply direct torque
        float torque = (Clockwise ? 1 : -1) * AngularForce * GetFalloff(body);
        body.ApplyTorque(torque * deltaTime);
    }
}

/// <summary>
/// Point effector - applies force from a single point (explosion)
/// </summary>
public class PointEffector : AreaEffector
{
    public float Magnitude { get; set; }
    public bool UseInverseSquare { get; set; } = true;

    public PointEffector(Vector2 position, float radius, float magnitude)
        : base(position, radius)
    {
        Magnitude = magnitude;
    }

    public override void ApplyEffect(RigidBody body, float deltaTime)
    {
        if (!IsBodyAffected(body))
            return;

        Vector2 direction = body.Position - Position;
        float distance = direction.Length;

        if (distance < 0.001f)
            return;

        direction = direction / distance;

        float forceMagnitude = Magnitude;

        if (UseInverseSquare)
        {
            forceMagnitude /= (distance * distance + 1.0f);
        }
        else
        {
            forceMagnitude *= GetFalloff(body);
        }

        body.ApplyForce(direction * forceMagnitude * deltaTime);
    }
}
