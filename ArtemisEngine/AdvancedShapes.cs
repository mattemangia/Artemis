namespace ArtemisEngine;

/// <summary>
/// Polygon shape for arbitrary convex polygons
/// Uses SAT for collision detection
/// </summary>
public class PolygonShape : Shape
{
    public Vector2[] LocalVertices { get; private set; }
    private Vector2[] _worldVertices;
    private Vector2[] _normals;

    public PolygonShape(Vector2[] vertices)
    {
        if (vertices.Length < 3)
            throw new ArgumentException("Polygon must have at least 3 vertices");

        LocalVertices = vertices;
        _worldVertices = new Vector2[vertices.Length];
        _normals = new Vector2[vertices.Length];

        // Pre-calculate normals
        for (int i = 0; i < vertices.Length; i++)
        {
            Vector2 edge = vertices[(i + 1) % vertices.Length] - vertices[i];
            _normals[i] = new Vector2(-edge.Y, edge.X).Normalized;
        }
    }

    public Vector2[] GetWorldVertices(Vector2 position, float rotation)
    {
        float cos = MathF.Cos(rotation);
        float sin = MathF.Sin(rotation);

        for (int i = 0; i < LocalVertices.Length; i++)
        {
            float x = LocalVertices[i].X * cos - LocalVertices[i].Y * sin;
            float y = LocalVertices[i].X * sin + LocalVertices[i].Y * cos;
            _worldVertices[i] = new Vector2(position.X + x, position.Y + y);
        }

        return _worldVertices;
    }

    public Vector2[] GetWorldNormals(float rotation)
    {
        float cos = MathF.Cos(rotation);
        float sin = MathF.Sin(rotation);

        for (int i = 0; i < _normals.Length; i++)
        {
            float x = _normals[i].X * cos - _normals[i].Y * sin;
            float y = _normals[i].X * sin + _normals[i].Y * cos;
            _normals[i] = new Vector2(x, y);
        }

        return _normals;
    }

    public override float CalculateInertia(float mass)
    {
        // Calculate moment of inertia for polygon
        float sum1 = 0;
        float sum2 = 0;

        for (int i = 0; i < LocalVertices.Length; i++)
        {
            Vector2 v1 = LocalVertices[i];
            Vector2 v2 = LocalVertices[(i + 1) % LocalVertices.Length];

            float cross = MathF.Abs(Vector2.Cross(v1, v2));
            sum1 += cross * (Vector2.Dot(v1, v1) + Vector2.Dot(v1, v2) + Vector2.Dot(v2, v2));
            sum2 += cross;
        }

        return (mass / 6.0f) * (sum1 / sum2);
    }

    /// <summary>
    /// Create a regular polygon (pentagon, hexagon, etc.)
    /// </summary>
    public static PolygonShape CreateRegular(int sides, float radius)
    {
        if (sides < 3)
            throw new ArgumentException("Must have at least 3 sides");

        Vector2[] vertices = new Vector2[sides];
        float angleStep = 2 * MathF.PI / sides;

        for (int i = 0; i < sides; i++)
        {
            float angle = i * angleStep;
            vertices[i] = new Vector2(
                radius * MathF.Cos(angle),
                radius * MathF.Sin(angle)
            );
        }

        return new PolygonShape(vertices);
    }
}

/// <summary>
/// Composite shape - combines multiple shapes into one body
/// Useful for complex objects
/// </summary>
public class CompoundShape : Shape
{
    public class ChildShape
    {
        public Shape Shape { get; set; }
        public Vector2 LocalPosition { get; set; }
        public float LocalRotation { get; set; }

        public ChildShape(Shape shape, Vector2 localPosition, float localRotation = 0)
        {
            Shape = shape;
            LocalPosition = localPosition;
            LocalRotation = localRotation;
        }
    }

    public List<ChildShape> Children { get; private set; }

    public CompoundShape()
    {
        Children = new List<ChildShape>();
    }

    public void AddChild(Shape shape, Vector2 localPosition, float localRotation = 0)
    {
        Children.Add(new ChildShape(shape, localPosition, localRotation));
    }

    public override float CalculateInertia(float mass)
    {
        // Distribute mass among children based on their "volume"
        float totalInertia = 0;
        float massPerChild = mass / Children.Count; // Simple distribution

        foreach (var child in Children)
        {
            float childInertia = child.Shape.CalculateInertia(massPerChild);

            // Parallel axis theorem: I = Icm + md^2
            float distance = child.LocalPosition.Length;
            totalInertia += childInertia + massPerChild * distance * distance;
        }

        return totalInertia;
    }

    /// <summary>
    /// Get AABB bounds for broadphase
    /// </summary>
    public (Vector2 min, Vector2 max) GetAABB(Vector2 position, float rotation)
    {
        Vector2 min = new Vector2(float.MaxValue, float.MaxValue);
        Vector2 max = new Vector2(float.MinValue, float.MinValue);

        foreach (var child in Children)
        {
            Vector2 worldPos = position + RotateVector(child.LocalPosition, rotation);
            float worldRot = rotation + child.LocalRotation;

            (Vector2 childMin, Vector2 childMax) = GetShapeAABB(child.Shape, worldPos, worldRot);

            min = new Vector2(Math.Min(min.X, childMin.X), Math.Min(min.Y, childMin.Y));
            max = new Vector2(Math.Max(max.X, childMax.X), Math.Max(max.Y, childMax.Y));
        }

        return (min, max);
    }

    private (Vector2, Vector2) GetShapeAABB(Shape shape, Vector2 position, float rotation)
    {
        if (shape is CircleShape circle)
        {
            return (
                new Vector2(position.X - circle.Radius, position.Y - circle.Radius),
                new Vector2(position.X + circle.Radius, position.Y + circle.Radius)
            );
        }
        else if (shape is BoxShape box)
        {
            // Simplified - use bounding circle for rotated box
            float diagonal = MathF.Sqrt(box.Width * box.Width + box.Height * box.Height) / 2;
            return (
                new Vector2(position.X - diagonal, position.Y - diagonal),
                new Vector2(position.X + diagonal, position.Y + diagonal)
            );
        }

        return (position, position);
    }

    private Vector2 RotateVector(Vector2 v, float angle)
    {
        float cos = MathF.Cos(angle);
        float sin = MathF.Sin(angle);
        return new Vector2(v.X * cos - v.Y * sin, v.X * sin + v.Y * cos);
    }
}

/// <summary>
/// Edge shape - for static terrain and platforms
/// Infinite length edge defined by two points
/// </summary>
public class EdgeShape : Shape
{
    public Vector2 Vertex1 { get; set; }
    public Vector2 Vertex2 { get; set; }
    public Vector2 Normal { get; private set; }

    public EdgeShape(Vector2 v1, Vector2 v2)
    {
        Vertex1 = v1;
        Vertex2 = v2;

        Vector2 edge = v2 - v1;
        Normal = new Vector2(-edge.Y, edge.X).Normalized;
    }

    public override float CalculateInertia(float mass)
    {
        // Edge is typically static, return 0
        return 0;
    }

    public Vector2 GetWorldNormal(float rotation)
    {
        float cos = MathF.Cos(rotation);
        float sin = MathF.Sin(rotation);

        return new Vector2(
            Normal.X * cos - Normal.Y * sin,
            Normal.X * sin + Normal.Y * cos
        );
    }
}
