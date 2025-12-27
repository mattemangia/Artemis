namespace ArtemisEngine;

/// <summary>
/// System for putting inactive bodies to sleep to improve performance
/// Sleeping bodies are not simulated until they receive an impulse or are near an awake body
/// </summary>
public class SleepingSystem
{
    private const float SleepVelocityThreshold = 0.1f;
    private const float SleepAngularVelocityThreshold = 0.05f;
    private const float TimeToSleep = 0.5f;

    public static void UpdateSleepState(RigidBody body, float deltaTime)
    {
        if (body.IsStatic || !body.CanSleep)
            return;

        float velocitySquared = body.Velocity.LengthSquared;
        float angularVelocitySquared = body.AngularVelocity * body.AngularVelocity;

        if (velocitySquared < SleepVelocityThreshold * SleepVelocityThreshold &&
            angularVelocitySquared < SleepAngularVelocityThreshold * SleepAngularVelocityThreshold)
        {
            body.SleepTimer += deltaTime;

            if (body.SleepTimer > TimeToSleep)
            {
                body.IsSleeping = true;
                body.Velocity = Vector2.Zero;
                body.AngularVelocity = 0;
            }
        }
        else
        {
            body.SleepTimer = 0;
            body.IsSleeping = false;
        }
    }

    public static void WakeUp(RigidBody body)
    {
        if (!body.IsStatic)
        {
            body.IsSleeping = false;
            body.SleepTimer = 0;
        }
    }

    public static void WakeUpNearby(PhysicsWorld world, RigidBody body, float radius = 5f)
    {
        Vector2 pos = body.Position;
        float radiusSquared = radius * radius;

        foreach (var other in world.Bodies)
        {
            if (other == body || other.IsStatic)
                continue;

            float distSquared = (other.Position - pos).LengthSquared;
            if (distSquared < radiusSquared)
            {
                WakeUp(other);
            }
        }
    }
}
