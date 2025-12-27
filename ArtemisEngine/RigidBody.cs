namespace ArtemisEngine;

public class RigidBody
{
    public Vector2 Position { get; set; }
    public Vector2 Velocity { get; set; }
    public float Rotation { get; set; }
    public float AngularVelocity { get; set; }

    public float Mass { get; set; }
    public float InverseMass { get; private set; }
    public float Restitution { get; set; } // Bounciness
    public float Friction { get; set; }

    public bool IsStatic { get; set; }
    public bool IsKinematic { get; set; }

    public Shape Shape { get; set; }

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
    }

    public void ApplyForce(Vector2 force)
    {
        if (IsStatic) return;
        Velocity += force * InverseMass;
    }

    public void ApplyImpulse(Vector2 impulse)
    {
        if (IsStatic) return;
        Velocity += impulse * InverseMass;
    }

    public void Update(float deltaTime)
    {
        if (IsStatic) return;

        Position += Velocity * deltaTime;
        Rotation += AngularVelocity * deltaTime;
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
