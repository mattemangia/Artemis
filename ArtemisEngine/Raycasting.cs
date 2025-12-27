namespace ArtemisEngine;

public struct Ray
{
    public Vector2 Origin { get; set; }
    public Vector2 Direction { get; set; }
    public float MaxDistance { get; set; }

    public Ray(Vector2 origin, Vector2 direction, float maxDistance = float.MaxValue)
    {
        Origin = origin;
        Direction = direction.Normalized;
        MaxDistance = maxDistance;
    }

    public Vector2 GetPoint(float distance)
    {
        return Origin + Direction * distance;
    }
}

public struct RaycastHit
{
    public RigidBody Body { get; set; }
    public Vector2 Point { get; set; }
    public Vector2 Normal { get; set; }
    public float Distance { get; set; }
    public bool Hit { get; set; }

    public static RaycastHit NoHit => new RaycastHit { Hit = false };
}

public static class Raycaster
{
    public static RaycastHit Raycast(Ray ray, IEnumerable<RigidBody> bodies, int layerMask = ~0)
    {
        RaycastHit closest = RaycastHit.NoHit;
        float closestDistance = ray.MaxDistance;

        foreach (var body in bodies)
        {
            // Skip if not in layer mask
            if ((body.CollisionLayer & layerMask) == 0)
                continue;

            RaycastHit hit;

            if (body.Shape is CircleShape circle)
            {
                hit = RaycastCircle(ray, body.Position, circle.Radius);
            }
            else if (body.Shape is BoxShape box)
            {
                hit = RaycastBox(ray, body.Position, box, body.Rotation);
            }
            else
            {
                continue;
            }

            if (hit.Hit && hit.Distance < closestDistance)
            {
                hit.Body = body;
                closest = hit;
                closestDistance = hit.Distance;
            }
        }

        return closest;
    }

    public static RaycastHit[] RaycastAll(Ray ray, IEnumerable<RigidBody> bodies, int layerMask = ~0)
    {
        List<RaycastHit> hits = new List<RaycastHit>();

        foreach (var body in bodies)
        {
            if ((body.CollisionLayer & layerMask) == 0)
                continue;

            RaycastHit hit;

            if (body.Shape is CircleShape circle)
            {
                hit = RaycastCircle(ray, body.Position, circle.Radius);
            }
            else if (body.Shape is BoxShape box)
            {
                hit = RaycastBox(ray, body.Position, box, body.Rotation);
            }
            else
            {
                continue;
            }

            if (hit.Hit && hit.Distance <= ray.MaxDistance)
            {
                hit.Body = body;
                hits.Add(hit);
            }
        }

        return hits.OrderBy(h => h.Distance).ToArray();
    }

    private static RaycastHit RaycastCircle(Ray ray, Vector2 center, float radius)
    {
        Vector2 m = ray.Origin - center;
        float b = Vector2.Dot(m, ray.Direction);
        float c = Vector2.Dot(m, m) - radius * radius;

        if (c > 0.0f && b > 0.0f)
            return RaycastHit.NoHit;

        float discriminant = b * b - c;

        if (discriminant < 0.0f)
            return RaycastHit.NoHit;

        float distance = -b - MathF.Sqrt(discriminant);

        if (distance < 0)
            distance = 0;

        Vector2 point = ray.Origin + ray.Direction * distance;
        Vector2 normal = (point - center).Normalized;

        return new RaycastHit
        {
            Hit = true,
            Point = point,
            Normal = normal,
            Distance = distance
        };
    }

    private static RaycastHit RaycastBox(Ray ray, Vector2 center, BoxShape box, float rotation)
    {
        // Transform ray to box local space
        float cos = MathF.Cos(-rotation);
        float sin = MathF.Sin(-rotation);

        Vector2 localOrigin = ray.Origin - center;
        Vector2 rotatedOrigin = new Vector2(
            localOrigin.X * cos - localOrigin.Y * sin,
            localOrigin.X * sin + localOrigin.Y * cos
        );

        Vector2 rotatedDirection = new Vector2(
            ray.Direction.X * cos - ray.Direction.Y * sin,
            ray.Direction.X * sin + ray.Direction.Y * cos
        );

        // AABB raycast in local space
        float halfW = box.Width / 2;
        float halfH = box.Height / 2;

        Vector2 min = new Vector2(-halfW, -halfH);
        Vector2 max = new Vector2(halfW, halfH);

        float tMin = 0;
        float tMax = ray.MaxDistance;

        for (int i = 0; i < 2; i++)
        {
            float origin = i == 0 ? rotatedOrigin.X : rotatedOrigin.Y;
            float dir = i == 0 ? rotatedDirection.X : rotatedDirection.Y;
            float minVal = i == 0 ? min.X : min.Y;
            float maxVal = i == 0 ? max.X : max.Y;

            if (MathF.Abs(dir) < 0.0001f)
            {
                if (origin < minVal || origin > maxVal)
                    return RaycastHit.NoHit;
            }
            else
            {
                float t1 = (minVal - origin) / dir;
                float t2 = (maxVal - origin) / dir;

                if (t1 > t2)
                {
                    float temp = t1;
                    t1 = t2;
                    t2 = temp;
                }

                tMin = Math.Max(tMin, t1);
                tMax = Math.Min(tMax, t2);

                if (tMin > tMax)
                    return RaycastHit.NoHit;
            }
        }

        Vector2 localPoint = rotatedOrigin + rotatedDirection * tMin;

        // Calculate normal in local space
        Vector2 localNormal;
        float epsilon = 0.0001f;

        if (MathF.Abs(localPoint.X - halfW) < epsilon)
            localNormal = new Vector2(1, 0);
        else if (MathF.Abs(localPoint.X + halfW) < epsilon)
            localNormal = new Vector2(-1, 0);
        else if (MathF.Abs(localPoint.Y - halfH) < epsilon)
            localNormal = new Vector2(0, 1);
        else if (MathF.Abs(localPoint.Y + halfH) < epsilon)
            localNormal = new Vector2(0, -1);
        else
            localNormal = new Vector2(0, 1);

        // Transform back to world space
        Vector2 worldNormal = new Vector2(
            localNormal.X * cos + localNormal.Y * sin,
            -localNormal.X * sin + localNormal.Y * cos
        );

        Vector2 worldPoint = ray.Origin + ray.Direction * tMin;

        return new RaycastHit
        {
            Hit = true,
            Point = worldPoint,
            Normal = worldNormal.Normalized,
            Distance = tMin
        };
    }
}
