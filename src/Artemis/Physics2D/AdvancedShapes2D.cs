using System;
using System.Collections.Generic;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Polygon shape for arbitrary convex polygons.
    /// Uses SAT (Separating Axis Theorem) for collision detection.
    /// </summary>
    public class PolygonShape2D : Shape2D
    {
        public Vector2D[] LocalVertices { get; private set; }
        private Vector2D[] _worldVertices;
        private Vector2D[] _normals;

        public PolygonShape2D(Vector2D[] vertices)
        {
            if (vertices.Length < 3)
                throw new ArgumentException("Polygon must have at least 3 vertices");

            LocalVertices = vertices;
            _worldVertices = new Vector2D[vertices.Length];
            _normals = new Vector2D[vertices.Length];

            // Pre-calculate normals
            for (int i = 0; i < vertices.Length; i++)
            {
                Vector2D edge = vertices[(i + 1) % vertices.Length] - vertices[i];
                _normals[i] = new Vector2D(-edge.Y, edge.X).Normalized;
            }
        }

        public Vector2D[] GetWorldVertices(Vector2D position, double rotation)
        {
            double cos = Math.Cos(rotation);
            double sin = Math.Sin(rotation);

            for (int i = 0; i < LocalVertices.Length; i++)
            {
                double x = LocalVertices[i].X * cos - LocalVertices[i].Y * sin;
                double y = LocalVertices[i].X * sin + LocalVertices[i].Y * cos;
                _worldVertices[i] = new Vector2D(position.X + x, position.Y + y);
            }

            return _worldVertices;
        }

        public Vector2D[] GetWorldNormals(double rotation)
        {
            double cos = Math.Cos(rotation);
            double sin = Math.Sin(rotation);

            var rotatedNormals = new Vector2D[_normals.Length];
            for (int i = 0; i < _normals.Length; i++)
            {
                double x = _normals[i].X * cos - _normals[i].Y * sin;
                double y = _normals[i].X * sin + _normals[i].Y * cos;
                rotatedNormals[i] = new Vector2D(x, y);
            }

            return rotatedNormals;
        }

        public override double CalculateInertia(double mass)
        {
            // Calculate moment of inertia for polygon
            double sum1 = 0;
            double sum2 = 0;

            for (int i = 0; i < LocalVertices.Length; i++)
            {
                Vector2D v1 = LocalVertices[i];
                Vector2D v2 = LocalVertices[(i + 1) % LocalVertices.Length];

                double cross = Math.Abs(Vector2D.Cross(v1, v2));
                sum1 += cross * (Vector2D.Dot(v1, v1) + Vector2D.Dot(v1, v2) + Vector2D.Dot(v2, v2));
                sum2 += cross;
            }

            return (mass / 6.0) * (sum1 / sum2);
        }

        public override double CalculateArea()
        {
            double area = 0;
            for (int i = 0; i < LocalVertices.Length; i++)
            {
                Vector2D v1 = LocalVertices[i];
                Vector2D v2 = LocalVertices[(i + 1) % LocalVertices.Length];
                area += Vector2D.Cross(v1, v2);
            }
            return Math.Abs(area) / 2.0;
        }

        /// <summary>
        /// Create a regular polygon (pentagon, hexagon, etc.).
        /// </summary>
        public static PolygonShape2D CreateRegular(int sides, double radius)
        {
            if (sides < 3)
                throw new ArgumentException("Must have at least 3 sides");

            Vector2D[] vertices = new Vector2D[sides];
            double angleStep = 2 * Math.PI / sides;

            for (int i = 0; i < sides; i++)
            {
                double angle = i * angleStep;
                vertices[i] = new Vector2D(
                    radius * Math.Cos(angle),
                    radius * Math.Sin(angle)
                );
            }

            return new PolygonShape2D(vertices);
        }
    }

    /// <summary>
    /// Composite shape - combines multiple shapes into one body.
    /// Useful for complex objects.
    /// </summary>
    public class CompoundShape2D : Shape2D
    {
        public class ChildShape
        {
            public Shape2D Shape { get; set; }
            public Vector2D LocalPosition { get; set; }
            public double LocalRotation { get; set; }

            public ChildShape(Shape2D shape, Vector2D localPosition, double localRotation = 0)
            {
                Shape = shape;
                LocalPosition = localPosition;
                LocalRotation = localRotation;
            }
        }

        public List<ChildShape> Children { get; private set; }

        public CompoundShape2D()
        {
            Children = new List<ChildShape>();
        }

        public void AddChild(Shape2D shape, Vector2D localPosition, double localRotation = 0)
        {
            Children.Add(new ChildShape(shape, localPosition, localRotation));
        }

        public override double CalculateInertia(double mass)
        {
            // Distribute mass among children based on their "volume"
            double totalInertia = 0;
            double massPerChild = mass / Children.Count; // Simple distribution

            foreach (var child in Children)
            {
                double childInertia = child.Shape.CalculateInertia(massPerChild);

                // Parallel axis theorem: I = Icm + md^2
                double distance = child.LocalPosition.Magnitude;
                totalInertia += childInertia + massPerChild * distance * distance;
            }

            return totalInertia;
        }

        public override double CalculateArea()
        {
            double totalArea = 0;
            foreach (var child in Children)
            {
                totalArea += child.Shape.CalculateArea();
            }
            return totalArea;
        }

        /// <summary>
        /// Get AABB bounds for broadphase.
        /// </summary>
        public (Vector2D min, Vector2D max) GetAABB(Vector2D position, double rotation)
        {
            Vector2D min = new Vector2D(double.MaxValue, double.MaxValue);
            Vector2D max = new Vector2D(double.MinValue, double.MinValue);

            foreach (var child in Children)
            {
                Vector2D worldPos = position + RotateVector(child.LocalPosition, rotation);
                double worldRot = rotation + child.LocalRotation;

                (Vector2D childMin, Vector2D childMax) = GetShapeAABB(child.Shape, worldPos, worldRot);

                min = new Vector2D(Math.Min(min.X, childMin.X), Math.Min(min.Y, childMin.Y));
                max = new Vector2D(Math.Max(max.X, childMax.X), Math.Max(max.Y, childMax.Y));
            }

            return (min, max);
        }

        private (Vector2D, Vector2D) GetShapeAABB(Shape2D shape, Vector2D position, double rotation)
        {
            if (shape is CircleShape circle)
            {
                return (
                    new Vector2D(position.X - circle.Radius, position.Y - circle.Radius),
                    new Vector2D(position.X + circle.Radius, position.Y + circle.Radius)
                );
            }
            else if (shape is BoxShape box)
            {
                // Simplified - use bounding circle for rotated box
                double diagonal = Math.Sqrt(box.Width * box.Width + box.Height * box.Height) / 2;
                return (
                    new Vector2D(position.X - diagonal, position.Y - diagonal),
                    new Vector2D(position.X + diagonal, position.Y + diagonal)
                );
            }

            return (position, position);
        }

        private Vector2D RotateVector(Vector2D v, double angle)
        {
            double cos = Math.Cos(angle);
            double sin = Math.Sin(angle);
            return new Vector2D(v.X * cos - v.Y * sin, v.X * sin + v.Y * cos);
        }
    }

    /// <summary>
    /// Edge shape - for static terrain and platforms.
    /// Infinite length edge defined by two points.
    /// </summary>
    public class EdgeShape2D : Shape2D
    {
        public Vector2D Vertex1 { get; set; }
        public Vector2D Vertex2 { get; set; }
        public Vector2D Normal { get; private set; }

        public EdgeShape2D(Vector2D v1, Vector2D v2)
        {
            Vertex1 = v1;
            Vertex2 = v2;

            Vector2D edge = v2 - v1;
            Normal = new Vector2D(-edge.Y, edge.X).Normalized;
        }

        public override double CalculateInertia(double mass)
        {
            // Edge is typically static, return 0
            return 0;
        }

        public override double CalculateArea()
        {
            // Edge has no area
            return 0;
        }

        public Vector2D GetWorldNormal(double rotation)
        {
            double cos = Math.Cos(rotation);
            double sin = Math.Sin(rotation);

            return new Vector2D(
                Normal.X * cos - Normal.Y * sin,
                Normal.X * sin + Normal.Y * cos
            );
        }

        public double GetLength()
        {
            return (Vertex2 - Vertex1).Magnitude;
        }
    }

    /// <summary>
    /// Chain shape - a series of connected edges forming a path.
    /// </summary>
    public class ChainShape2D : Shape2D
    {
        public Vector2D[] Vertices { get; private set; }
        public bool IsLoop { get; private set; }

        public ChainShape2D(Vector2D[] vertices, bool isLoop = false)
        {
            if (vertices.Length < 2)
                throw new ArgumentException("Chain must have at least 2 vertices");

            Vertices = vertices;
            IsLoop = isLoop;
        }

        public EdgeShape2D GetEdge(int index)
        {
            int next = (index + 1) % Vertices.Length;
            return new EdgeShape2D(Vertices[index], Vertices[next]);
        }

        public int EdgeCount => IsLoop ? Vertices.Length : Vertices.Length - 1;

        public override double CalculateInertia(double mass)
        {
            return 0; // Chains are typically static
        }

        public override double CalculateArea()
        {
            return 0; // Chain has no area
        }
    }
}
