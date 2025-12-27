namespace ArtemisEngine;

public class RigidBody
{
    public Vector2 Position { get; set; }
    public Vector2 Velocity { get; set; }
    public float Rotation { get; set; }
    public float AngularVelocity { get; set; }

    public float Mass { get; set; }
    public float InverseMass { get; private set; }
    public float Inertia { get; private set; }
    public float InverseInertia { get; private set; }

    public float Restitution { get; set; } // Bounciness
    public float Friction { get; set; }

    public bool IsStatic { get; set; }
    public bool IsKinematic { get; set; }
    public bool IsTrigger { get; set; } // Sensor - no physical collision

    // Sleeping system
    public bool IsSleeping { get; set; }
    public bool CanSleep { get; set; } = true;
    public float SleepTimer { get; set; }

    // Collision filtering
    public int CollisionLayer { get; set; } = CollisionLayers.Default;
    public int CollisionMask { get; set; } = CollisionLayers.Everything;

    // Collision callbacks
    public ICollisionListener? CollisionListener { get; set; }

    public Shape Shape { get; set; }

    // User data for game-specific information
    public object? UserData { get; set; }

    public RigidBody(Vector2 position, float mass, Shape shape, bool isStatic = false)
    {
        Position = position;
        Mass = mass;
        InverseMass = isStatic ? 0 : (mass > 0 ? 1.0f / mass : 0);
        Shape = shape;
        IsStatic = isStatic;

        Velocity = Vector2.Zero;
        AngularVelocity = 0;
        Rotation = 0;
        Restitution = 0.5f;
        Friction = 0.3f;

        // Calculate inertia
        Inertia = shape.CalculateInertia(mass);
        InverseInertia = Inertia > 0 ? 1.0f / Inertia : 0;

        IsSleeping = false;
        SleepTimer = 0;
    }

    public void ApplyForce(Vector2 force)
    {
        if (IsStatic || IsSleeping) return;
        Velocity += force * InverseMass;
        WakeUp();
    }

    public void ApplyImpulse(Vector2 impulse)
    {
        if (IsStatic || IsSleeping) return;
        Velocity += impulse * InverseMass;
        WakeUp();
    }

    public void ApplyImpulseAtPoint(Vector2 impulse, Vector2 point)
    {
        if (IsStatic || IsSleeping) return;

        Velocity += impulse * InverseMass;

        Vector2 r = point - Position;
        float torque = Vector2.Cross(r, impulse);
        AngularVelocity += torque * InverseInertia;

        WakeUp();
    }

    public void ApplyTorque(float torque)
    {
        if (IsStatic || IsSleeping) return;
        AngularVelocity += torque * InverseInertia;
        WakeUp();
    }

    public void WakeUp()
    {
        if (!IsStatic)
        {
            IsSleeping = false;
            SleepTimer = 0;
        }
    }

    public void Update(float deltaTime)
    {
        if (IsStatic || IsSleeping) return;

        Position += Velocity * deltaTime;
        Rotation += AngularVelocity * deltaTime;
    }

    public bool CanCollideWith(RigidBody other)
    {
        // Static bodies don't collide with each other
        if (IsStatic && other.IsStatic)
            return false;

        // Check layer filtering
        if ((CollisionLayer & other.CollisionMask) == 0)
            return false;
        if ((other.CollisionLayer & CollisionMask) == 0)
            return false;

        return true;
    }
}

public abstract class Shape
{
    public abstract float CalculateInertia(float mass);
}

public class CircleShape : Shape
{
    public float Radius { get; set; }

    public CircleShape(float radius)
    {
        Radius = radius;
    }

    public override float CalculateInertia(float mass)
    {
        return 0.5f * mass * Radius * Radius;
    }
}

public class BoxShape : Shape
{
    public float Width { get; set; }
    public float Height { get; set; }

    public BoxShape(float width, float height)
    {
        Width = width;
        Height = height;
    }

    public override float CalculateInertia(float mass)
    {
        return mass * (Width * Width + Height * Height) / 12.0f;
    }

    public Vector2[] GetVertices(Vector2 position, float rotation)
    {
        float halfW = Width / 2;
        float halfH = Height / 2;

        Vector2[] localVertices = new[]
        {
            new Vector2(-halfW, -halfH),
            new Vector2(halfW, -halfH),
            new Vector2(halfW, halfH),
            new Vector2(-halfW, halfH)
        };

        Vector2[] worldVertices = new Vector2[4];
        float cos = MathF.Cos(rotation);
        float sin = MathF.Sin(rotation);

        for (int i = 0; i < 4; i++)
        {
            float x = localVertices[i].X * cos - localVertices[i].Y * sin;
            float y = localVertices[i].X * sin + localVertices[i].Y * cos;
            worldVertices[i] = new Vector2(position.X + x, position.Y + y);
        }

        return worldVertices;
    }
}
