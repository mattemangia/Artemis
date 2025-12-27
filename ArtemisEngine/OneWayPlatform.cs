namespace ArtemisEngine;

/// <summary>
/// One-way platform collision - objects can pass through from one direction
/// Commonly used for platforms you can jump up through
/// </summary>
public class OneWayPlatformBehavior
{
    public Vector2 AllowedDirection { get; set; } // Direction you CAN pass through
    public float Threshold { get; set; } = 0.1f; // Velocity threshold

    public OneWayPlatformBehavior(Vector2 allowedDirection)
    {
        AllowedDirection = allowedDirection.Normalized;
    }

    /// <summary>
    /// Check if collision should be processed based on approach direction
    /// </summary>
    public bool ShouldCollide(RigidBody moving, RigidBody platform, Vector2 collisionNormal)
    {
        // If moving in the allowed direction, pass through
        float approachDot = Vector2.Dot(moving.Velocity.Normalized, AllowedDirection);

        if (approachDot > Threshold)
        {
            // Moving in allowed direction - no collision
            return false;
        }

        // Check if approaching from the correct side
        float normalDot = Vector2.Dot(collisionNormal, AllowedDirection);

        if (normalDot > 0)
        {
            // Collision normal points in allowed direction - allow pass through
            return false;
        }

        // Otherwise, collide normally
        return true;
    }

    /// <summary>
    /// Create a standard "jump-through" platform (pass through from below)
    /// </summary>
    public static OneWayPlatformBehavior CreateJumpThrough()
    {
        return new OneWayPlatformBehavior(new Vector2(0, 1)); // Allow upward passage
    }

    /// <summary>
    /// Create a fall-through platform (pass through from above)
    /// </summary>
    public static OneWayPlatformBehavior CreateFallThrough()
    {
        return new OneWayPlatformBehavior(new Vector2(0, -1)); // Allow downward passage
    }
}

/// <summary>
/// Extended RigidBody properties for one-way platforms
/// </summary>
public static class OneWayPlatformExtensions
{
    private static Dictionary<RigidBody, OneWayPlatformBehavior> _platformBehaviors = new();

    public static void SetOneWayPlatform(this RigidBody body, OneWayPlatformBehavior behavior)
    {
        _platformBehaviors[body] = behavior;
    }

    public static OneWayPlatformBehavior? GetOneWayPlatform(this RigidBody body)
    {
        return _platformBehaviors.TryGetValue(body, out var behavior) ? behavior : null;
    }

    public static bool IsOneWayPlatform(this RigidBody body)
    {
        return _platformBehaviors.ContainsKey(body);
    }

    public static void RemoveOneWayPlatform(this RigidBody body)
    {
        _platformBehaviors.Remove(body);
    }

    /// <summary>
    /// Check if collision between two bodies should be ignored due to one-way platform
    /// </summary>
    public static bool ShouldIgnoreCollision(RigidBody bodyA, RigidBody bodyB, Vector2 normal)
    {
        var platformA = bodyA.GetOneWayPlatform();
        var platformB = bodyB.GetOneWayPlatform();

        if (platformA != null && !platformA.ShouldCollide(bodyB, bodyA, normal))
            return true;

        if (platformB != null && !platformB.ShouldCollide(bodyA, bodyB, -normal))
            return true;

        return false;
    }
}
