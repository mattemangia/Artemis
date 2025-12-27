namespace ArtemisEngine;

/// <summary>
/// Advanced collision detection algorithms
/// </summary>
public static class CollisionDetection
{
    /// <summary>
    /// SAT (Separating Axis Theorem) collision detection for oriented boxes
    /// Supports full rotation - no more AABB limitations!
    /// </summary>
    public static bool BoxVsBoxSAT(RigidBody bodyA, RigidBody bodyB, BoxShape boxA, BoxShape boxB, out Collision collision)
    {
        collision = new Collision();

        // Get vertices in world space
        Vector2[] verticesA = boxA.GetVertices(bodyA.Position, bodyA.Rotation);
        Vector2[] verticesB = boxB.GetVertices(bodyB.Position, bodyB.Rotation);

        // Test all axes (4 edge normals total: 2 from each box)
        Vector2[] axes = new Vector2[4];

        // Axes from box A
        axes[0] = GetEdgeNormal(verticesA[0], verticesA[1]);
        axes[1] = GetEdgeNormal(verticesA[1], verticesA[2]);

        // Axes from box B
        axes[2] = GetEdgeNormal(verticesB[0], verticesB[1]);
        axes[3] = GetEdgeNormal(verticesB[1], verticesB[2]);

        float minOverlap = float.MaxValue;
        Vector2 smallestAxis = Vector2.Zero;

        // Test each axis
        foreach (var axis in axes)
        {
            // Project both shapes onto the axis
            ProjectVertices(verticesA, axis, out float minA, out float maxA);
            ProjectVertices(verticesB, axis, out float minB, out float maxB);

            // Check for separation
            if (maxA < minB || maxB < minA)
            {
                // Shapes are separated on this axis
                return false;
            }

            // Calculate overlap
            float overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);

            if (overlap < minOverlap)
            {
                minOverlap = overlap;
                smallestAxis = axis;
            }
        }

        // If we get here, there's a collision
        collision.Penetration = minOverlap;

        // Ensure normal points from A to B
        Vector2 direction = bodyB.Position - bodyA.Position;
        if (Vector2.Dot(direction, smallestAxis) < 0)
        {
            smallestAxis = -smallestAxis;
        }

        collision.Normal = smallestAxis;

        // Calculate contact point (average of penetrating vertices)
        collision.ContactPoint = CalculateContactPoint(verticesA, verticesB, smallestAxis);

        return true;
    }

    /// <summary>
    /// Improved Circle vs Box collision with proper rotation support
    /// </summary>
    public static bool CircleVsBoxImproved(RigidBody circleBody, RigidBody boxBody, CircleShape circle, BoxShape box, out Collision collision)
    {
        collision = new Collision();

        // Transform circle center to box local space
        Vector2 localCirclePos = WorldToLocal(circleBody.Position, boxBody.Position, boxBody.Rotation);

        // Clamp to box bounds
        float halfW = box.Width / 2;
        float halfH = box.Height / 2;

        Vector2 closest = new Vector2(
            Math.Clamp(localCirclePos.X, -halfW, halfW),
            Math.Clamp(localCirclePos.Y, -halfH, halfH)
        );

        Vector2 localPoint = localCirclePos - closest;
        float distanceSquared = localPoint.LengthSquared;

        if (distanceSquared >= circle.Radius * circle.Radius)
            return false;

        float distance = MathF.Sqrt(distanceSquared);

        // Calculate normal in local space
        Vector2 localNormal;
        if (distance > 0.0001f)
        {
            localNormal = localPoint / distance;
        }
        else
        {
            // Circle center inside box - find closest edge
            float[] distances = {
                halfW - MathF.Abs(localCirclePos.X),
                halfH - MathF.Abs(localCirclePos.Y)
            };

            if (distances[0] < distances[1])
            {
                localNormal = new Vector2(localCirclePos.X > 0 ? 1 : -1, 0);
            }
            else
            {
                localNormal = new Vector2(0, localCirclePos.Y > 0 ? 1 : -1);
            }
        }

        // Transform normal back to world space
        collision.Normal = LocalToWorldDirection(localNormal, boxBody.Rotation);
        collision.Penetration = circle.Radius - distance;

        // Contact point in world space
        Vector2 localContactPoint = closest;
        collision.ContactPoint = LocalToWorld(localContactPoint, boxBody.Position, boxBody.Rotation);

        return true;
    }

    /// <summary>
    /// Calculate accurate contact point between two boxes
    /// </summary>
    private static Vector2 CalculateContactPoint(Vector2[] verticesA, Vector2[] verticesB, Vector2 axis)
    {
        // Find vertices on the penetration side
        List<Vector2> contactPoints = new List<Vector2>();

        // Check which vertices of A are inside B
        foreach (var vertex in verticesA)
        {
            if (IsPointInPolygon(vertex, verticesB))
            {
                contactPoints.Add(vertex);
            }
        }

        // Check which vertices of B are inside A
        foreach (var vertex in verticesB)
        {
            if (IsPointInPolygon(vertex, verticesA))
            {
                contactPoints.Add(vertex);
            }
        }

        // If we have contact points, average them
        if (contactPoints.Count > 0)
        {
            Vector2 sum = Vector2.Zero;
            foreach (var point in contactPoints)
            {
                sum += point;
            }
            return sum / contactPoints.Count;
        }

        // Fallback: midpoint between centers
        Vector2 centerA = (verticesA[0] + verticesA[2]) / 2;
        Vector2 centerB = (verticesB[0] + verticesB[2]) / 2;
        return (centerA + centerB) / 2;
    }

    private static bool IsPointInPolygon(Vector2 point, Vector2[] polygon)
    {
        bool inside = false;
        int j = polygon.Length - 1;

        for (int i = 0; i < polygon.Length; j = i++)
        {
            if (((polygon[i].Y > point.Y) != (polygon[j].Y > point.Y)) &&
                (point.X < (polygon[j].X - polygon[i].X) * (point.Y - polygon[i].Y) /
                 (polygon[j].Y - polygon[i].Y) + polygon[i].X))
            {
                inside = !inside;
            }
        }

        return inside;
    }

    private static Vector2 GetEdgeNormal(Vector2 p1, Vector2 p2)
    {
        Vector2 edge = p2 - p1;
        Vector2 normal = new Vector2(-edge.Y, edge.X);
        return normal.Normalized;
    }

    private static void ProjectVertices(Vector2[] vertices, Vector2 axis, out float min, out float max)
    {
        min = float.MaxValue;
        max = float.MinValue;

        foreach (var vertex in vertices)
        {
            float projection = Vector2.Dot(vertex, axis);
            min = Math.Min(min, projection);
            max = Math.Max(max, projection);
        }
    }

    private static Vector2 WorldToLocal(Vector2 worldPos, Vector2 origin, float rotation)
    {
        Vector2 delta = worldPos - origin;
        float cos = MathF.Cos(-rotation);
        float sin = MathF.Sin(-rotation);

        return new Vector2(
            delta.X * cos - delta.Y * sin,
            delta.X * sin + delta.Y * cos
        );
    }

    private static Vector2 LocalToWorld(Vector2 localPos, Vector2 origin, float rotation)
    {
        float cos = MathF.Cos(rotation);
        float sin = MathF.Sin(rotation);

        return new Vector2(
            origin.X + localPos.X * cos - localPos.Y * sin,
            origin.Y + localPos.X * sin + localPos.Y * cos
        );
    }

    private static Vector2 LocalToWorldDirection(Vector2 localDir, float rotation)
    {
        float cos = MathF.Cos(rotation);
        float sin = MathF.Sin(rotation);

        return new Vector2(
            localDir.X * cos - localDir.Y * sin,
            localDir.X * sin + localDir.Y * cos
        );
    }
}
