namespace ArtemisEngine;

/// <summary>
/// Continuous Collision Detection (CCD) to prevent tunneling for fast-moving objects
/// Uses swept collision detection and time-of-impact calculation
/// </summary>
public static class ContinuousCollision
{
    /// <summary>
    /// Perform swept collision detection for a moving body
    /// Returns true if collision occurs within the timestep
    /// </summary>
    public static bool SweptCollision(RigidBody moving, RigidBody target, float deltaTime, out float timeOfImpact, out Collision collision)
    {
        timeOfImpact = 1.0f;
        collision = new Collision();

        // Calculate trajectory
        Vector2 startPos = moving.Position;
        Vector2 endPos = moving.Position + moving.Velocity * deltaTime;
        Vector2 trajectory = endPos - startPos;

        if (trajectory.LengthSquared < 0.0001f)
        {
            // Not moving significantly, use discrete collision
            return false;
        }

        // Use conservative advancement or raycast-based CCD depending on shape
        if (moving.Shape is CircleShape movingCircle && target.Shape is CircleShape targetCircle)
        {
            return SweptCircleVsCircle(startPos, endPos, movingCircle, target.Position, targetCircle, out timeOfImpact, out collision);
        }
        else if (moving.Shape is CircleShape circle)
        {
            // Circle vs anything - use raycast
            return RaycastCCD(startPos, trajectory, circle.Radius, target, out timeOfImpact, out collision);
        }
        else if (moving.Shape is BoxShape box)
        {
            // Box - use conservative advancement with multiple samples
            return ConservativeAdvancement(moving, target, deltaTime, out timeOfImpact, out collision);
        }

        return false;
    }

    private static bool SweptCircleVsCircle(Vector2 startA, Vector2 endA, CircleShape circleA,
                                           Vector2 posB, CircleShape circleB,
                                           out float toi, out Collision collision)
    {
        toi = 1.0f;
        collision = new Collision();

        Vector2 trajectory = endA - startA;
        Vector2 centerToCenter = posB - startA;
        float combinedRadius = circleA.Radius + circleB.Radius;

        // Solve quadratic equation for time of impact
        // |startA + t*trajectory - posB|^2 = combinedRadius^2
        // Let d = startA - posB = -centerToCenter
        // |d + t*trajectory|^2 = combinedRadius^2
        // Expanding: t^2*|trajectory|^2 + 2*t*(d·trajectory) + |d|^2 = combinedRadius^2
        // So: a*t^2 + b*t + c = 0 where b = 2*(d·trajectory) = -2*(centerToCenter·trajectory)

        float a = trajectory.LengthSquared;
        float b = -2 * Vector2.Dot(trajectory, centerToCenter);
        float c = centerToCenter.LengthSquared - combinedRadius * combinedRadius;

        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
        {
            // No collision
            return false;
        }

        float sqrtDisc = MathF.Sqrt(discriminant);
        float t1 = (-b - sqrtDisc) / (2 * a);
        float t2 = (-b + sqrtDisc) / (2 * a);

        // We want the first impact (smallest positive t)
        if (t1 >= 0 && t1 <= 1.0f)
        {
            toi = t1;
        }
        else if (t2 >= 0 && t2 <= 1.0f)
        {
            toi = t2;
        }
        else
        {
            return false;
        }

        // Calculate collision info at time of impact
        Vector2 impactPosA = startA + trajectory * toi;
        Vector2 delta = posB - impactPosA;
        float distance = delta.Length;

        collision.Normal = distance > 0 ? delta / distance : new Vector2(1, 0);
        collision.Penetration = 0; // No penetration at TOI
        collision.ContactPoint = impactPosA + collision.Normal * circleA.Radius;

        return true;
    }

    private static bool RaycastCCD(Vector2 start, Vector2 trajectory, float radius, RigidBody target,
                                  out float toi, out Collision collision)
    {
        toi = 1.0f;
        collision = new Collision();

        // Raycast along the trajectory
        Ray ray = new Ray(start, trajectory.Normalized, trajectory.Length + radius);

        // Create a temporary physics world just for raycasting (hacky but works)
        // In production, you'd pass the world as parameter
        var hit = Raycaster.Raycast(ray, new[] { target }, ~0);

        if (hit.Hit && hit.Distance < trajectory.Length)
        {
            toi = hit.Distance / trajectory.Length;
            collision.Normal = hit.Normal;
            collision.ContactPoint = hit.Point;
            collision.Penetration = 0;
            return true;
        }

        return false;
    }

    private static bool ConservativeAdvancement(RigidBody moving, RigidBody target, float deltaTime,
                                               out float toi, out Collision collision)
    {
        toi = 1.0f;
        collision = new Collision();

        const int maxSteps = 10;
        const float tolerance = 0.01f;

        Vector2 startPos = moving.Position;
        Vector2 trajectory = moving.Velocity * deltaTime;
        float currentT = 0;

        // Incrementally advance until collision or end of trajectory
        for (int step = 0; step < maxSteps; step++)
        {
            float remainingT = 1.0f - currentT;
            if (remainingT <= 0)
                break;

            // Test current position
            moving.Position = startPos + trajectory * currentT;

            bool collides = TestCollision(moving, target, out var testCollision);

            if (collides && testCollision.Penetration > tolerance)
            {
                // Found collision, back up to find exact TOI
                toi = currentT;
                collision = testCollision;
                moving.Position = startPos; // Restore position
                return true;
            }

            // Advance conservatively
            float stepSize = collides ? tolerance / 2 : remainingT / 2;
            currentT += stepSize;
        }

        moving.Position = startPos; // Restore position
        return false;
    }

    private static bool TestCollision(RigidBody bodyA, RigidBody bodyB, out Collision collision)
    {
        collision = new Collision();

        if (bodyA.Shape is CircleShape circleA && bodyB.Shape is CircleShape circleB)
        {
            Vector2 delta = bodyB.Position - bodyA.Position;
            float distSq = delta.LengthSquared;
            float radiusSum = circleA.Radius + circleB.Radius;

            if (distSq < radiusSum * radiusSum)
            {
                float dist = MathF.Sqrt(distSq);
                collision.Normal = dist > 0 ? delta / dist : new Vector2(1, 0);
                collision.Penetration = radiusSum - dist;
                return true;
            }
        }
        else if (bodyA.Shape is BoxShape boxA && bodyB.Shape is BoxShape boxB)
        {
            return CollisionDetection.BoxVsBoxSAT(bodyA, bodyB, boxA, boxB, out collision);
        }
        else if (bodyA.Shape is CircleShape circle && bodyB.Shape is BoxShape box)
        {
            return CollisionDetection.CircleVsBoxImproved(bodyA, bodyB, circle, box, out collision);
        }
        else if (bodyA.Shape is BoxShape box2 && bodyB.Shape is CircleShape circle2)
        {
            bool result = CollisionDetection.CircleVsBoxImproved(bodyB, bodyA, circle2, box2, out collision);
            if (result)
                collision.Normal = -collision.Normal;
            return result;
        }

        return false;
    }

    /// <summary>
    /// Determine if a body should use CCD based on its velocity
    /// </summary>
    public static bool ShouldUseCCD(RigidBody body, float deltaTime)
    {
        if (body.IsStatic || !body.UseCCD)
            return false;

        // Use CCD if object moves more than its smallest dimension in one frame
        float movement = body.Velocity.Length * deltaTime;
        float minSize = GetMinSize(body.Shape);

        return movement > minSize * 0.5f; // Threshold: half the body size
    }

    private static float GetMinSize(Shape shape)
    {
        if (shape is CircleShape circle)
            return circle.Radius * 2;
        else if (shape is BoxShape box)
            return Math.Min(box.Width, box.Height);

        return 1.0f;
    }
}
