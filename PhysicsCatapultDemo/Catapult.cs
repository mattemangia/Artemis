using ArtemisEngine;

namespace PhysicsCatapultDemo;

public class Catapult
{
    public Vector2 Position { get; set; }
    public Vector2 AimDirection { get; set; }
    public float Power { get; set; }
    public float MinPower { get; } = 100f;
    public float MaxPower { get; } = 800f;

    public Catapult(Vector2 position)
    {
        Position = position;
        AimDirection = new Vector2(1, -1).Normalized;
        Power = 400f;
    }

    public void SetAim(Vector2 targetPosition)
    {
        Vector2 direction = targetPosition - Position;
        if (direction.LengthSquared > 0.01f)
        {
            AimDirection = direction.Normalized;
        }
    }

    public void SetAimAngle(float angleDegrees)
    {
        float angleRadians = angleDegrees * MathF.PI / 180f;
        AimDirection = new Vector2(MathF.Cos(angleRadians), -MathF.Sin(angleRadians));
    }

    public void AdjustPower(float delta)
    {
        Power = Math.Clamp(Power + delta, MinPower, MaxPower);
    }

    public GameObject LaunchProjectile(float radius = 0.5f, float mass = 1.0f)
    {
        Vector2 spawnPos = Position + AimDirection * 2f;
        Vector2 velocity = AimDirection * Power;

        var shape = new CircleShape(radius);
        var body = new RigidBody(spawnPos, mass, shape)
        {
            Velocity = velocity,
            Restitution = 0.6f,
            Friction = 0.3f
        };

        return new GameObject(body, GameObjectType.Projectile);
    }

    public Vector2 GetAimEndPoint(float length = 5f)
    {
        return Position + AimDirection * length;
    }
}
