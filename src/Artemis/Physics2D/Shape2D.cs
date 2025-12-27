using System;
using System.Runtime.CompilerServices;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Base class for 2D collision shapes.
    /// </summary>
    public abstract class Shape2D
    {
        /// <summary>
        /// Gets the type of this shape.
        /// </summary>
        public abstract ShapeType2D Type { get; }

        /// <summary>
        /// Computes the axis-aligned bounding box for this shape.
        /// </summary>
        public abstract AABB2D ComputeAABB(Vector2D position, double rotation);

        /// <summary>
        /// Computes the mass properties for this shape given a density.
        /// </summary>
        public abstract MassData2D ComputeMass(double density);
    }

    /// <summary>
    /// Defines 2D shape types.
    /// </summary>
    public enum ShapeType2D
    {
        Circle,
        Box,
        Polygon,
        Capsule,
        Edge,
        Chain
    }

    /// <summary>
    /// Mass and inertia data for a 2D shape.
    /// </summary>
    public struct MassData2D
    {
        /// <summary>Mass in kg.</summary>
        public double Mass;

        /// <summary>Moment of inertia around the center of mass.</summary>
        public double Inertia;

        /// <summary>Center of mass in local coordinates.</summary>
        public Vector2D Center;
    }

    /// <summary>
    /// 2D axis-aligned bounding box.
    /// </summary>
    public struct AABB2D
    {
        public Vector2D Min;
        public Vector2D Max;

        public AABB2D(Vector2D min, Vector2D max)
        {
            Min = min;
            Max = max;
        }

        public readonly Vector2D Center => (Min + Max) * 0.5;
        public readonly Vector2D Extents => (Max - Min) * 0.5;
        public readonly double Width => Max.X - Min.X;
        public readonly double Height => Max.Y - Min.Y;
        public readonly double Area => Width * Height;
        public readonly double Perimeter => 2 * (Width + Height);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Contains(Vector2D point)
            => point.X >= Min.X && point.X <= Max.X &&
               point.Y >= Min.Y && point.Y <= Max.Y;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Overlaps(AABB2D other)
            => Min.X <= other.Max.X && Max.X >= other.Min.X &&
               Min.Y <= other.Max.Y && Max.Y >= other.Min.Y;

        public static AABB2D Combine(AABB2D a, AABB2D b)
            => new(Vector2D.Min(a.Min, b.Min), Vector2D.Max(a.Max, b.Max));

        public void Expand(double amount)
        {
            var delta = new Vector2D(amount);
            Min -= delta;
            Max += delta;
        }
    }

    /// <summary>
    /// Circle shape for 2D physics.
    /// </summary>
    public class CircleShape : Shape2D
    {
        /// <summary>Radius of the circle.</summary>
        public double Radius { get; set; } = 0.5;

        /// <summary>Local offset from body center.</summary>
        public Vector2D Offset { get; set; } = Vector2D.Zero;

        public override ShapeType2D Type => ShapeType2D.Circle;

        public CircleShape() { }

        public CircleShape(double radius, Vector2D? offset = null)
        {
            Radius = radius;
            Offset = offset ?? Vector2D.Zero;
        }

        public override AABB2D ComputeAABB(Vector2D position, double rotation)
        {
            var center = position + Vector2D.Rotate(Offset, rotation);
            var r = new Vector2D(Radius);
            return new AABB2D(center - r, center + r);
        }

        public override MassData2D ComputeMass(double density)
        {
            double mass = density * Math.PI * Radius * Radius;
            // I = (1/2) * m * r^2
            double inertia = 0.5 * mass * Radius * Radius;

            // Add parallel axis theorem for offset
            if (Offset.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                inertia += mass * Offset.MagnitudeSquared;
            }

            return new MassData2D
            {
                Mass = mass,
                Inertia = inertia,
                Center = Offset
            };
        }
    }

    /// <summary>
    /// Box (rectangle) shape for 2D physics.
    /// </summary>
    public class BoxShape : Shape2D
    {
        /// <summary>Half-width of the box.</summary>
        public double HalfWidth { get; set; } = 0.5;

        /// <summary>Half-height of the box.</summary>
        public double HalfHeight { get; set; } = 0.5;

        /// <summary>Local offset from body center.</summary>
        public Vector2D Offset { get; set; } = Vector2D.Zero;

        public double Width => HalfWidth * 2;
        public double Height => HalfHeight * 2;

        public override ShapeType2D Type => ShapeType2D.Box;

        public BoxShape() { }

        public BoxShape(double halfWidth, double halfHeight, Vector2D? offset = null)
        {
            HalfWidth = halfWidth;
            HalfHeight = halfHeight;
            Offset = offset ?? Vector2D.Zero;
        }

        public BoxShape(double size) : this(size * 0.5, size * 0.5) { }

        /// <summary>
        /// Gets the four corners of the box in local space.
        /// </summary>
        public Vector2D[] GetVertices()
        {
            return new[]
            {
                Offset + new Vector2D(-HalfWidth, -HalfHeight),
                Offset + new Vector2D(HalfWidth, -HalfHeight),
                Offset + new Vector2D(HalfWidth, HalfHeight),
                Offset + new Vector2D(-HalfWidth, HalfHeight)
            };
        }

        /// <summary>
        /// Gets the four corners of the box in world space.
        /// </summary>
        public Vector2D[] GetWorldVertices(Vector2D position, double rotation)
        {
            var vertices = GetVertices();
            for (int i = 0; i < vertices.Length; i++)
            {
                vertices[i] = position + Vector2D.Rotate(vertices[i], rotation);
            }
            return vertices;
        }

        public override AABB2D ComputeAABB(Vector2D position, double rotation)
        {
            var vertices = GetWorldVertices(position, rotation);
            var min = vertices[0];
            var max = vertices[0];

            for (int i = 1; i < vertices.Length; i++)
            {
                min = Vector2D.Min(min, vertices[i]);
                max = Vector2D.Max(max, vertices[i]);
            }

            return new AABB2D(min, max);
        }

        public override MassData2D ComputeMass(double density)
        {
            double mass = density * Width * Height;
            // I = (1/12) * m * (w^2 + h^2)
            double inertia = mass * (Width * Width + Height * Height) / 12.0;

            // Add parallel axis theorem for offset
            if (Offset.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                inertia += mass * Offset.MagnitudeSquared;
            }

            return new MassData2D
            {
                Mass = mass,
                Inertia = inertia,
                Center = Offset
            };
        }
    }

    /// <summary>
    /// Convex polygon shape for 2D physics.
    /// </summary>
    public class PolygonShape : Shape2D
    {
        private Vector2D[] _vertices;
        private Vector2D[] _normals;

        /// <summary>Gets the vertices of the polygon in local space.</summary>
        public Vector2D[] Vertices => _vertices;

        /// <summary>Gets the edge normals of the polygon.</summary>
        public Vector2D[] Normals => _normals;

        /// <summary>Number of vertices.</summary>
        public int VertexCount => _vertices?.Length ?? 0;

        public override ShapeType2D Type => ShapeType2D.Polygon;

        public PolygonShape()
        {
            _vertices = Array.Empty<Vector2D>();
            _normals = Array.Empty<Vector2D>();
        }

        /// <summary>
        /// Creates a polygon from vertices (must be convex, counter-clockwise winding).
        /// </summary>
        public PolygonShape(Vector2D[] vertices)
        {
            _vertices = Array.Empty<Vector2D>();
            _normals = Array.Empty<Vector2D>();
            SetVertices(vertices);
        }

        /// <summary>
        /// Sets the vertices of the polygon.
        /// </summary>
        public void SetVertices(Vector2D[] vertices)
        {
            if (vertices.Length < 3)
                throw new ArgumentException("Polygon must have at least 3 vertices");

            _vertices = new Vector2D[vertices.Length];
            _normals = new Vector2D[vertices.Length];

            Array.Copy(vertices, _vertices, vertices.Length);
            ComputeNormals();
        }

        /// <summary>
        /// Creates a regular polygon with the given number of sides and radius.
        /// </summary>
        public static PolygonShape CreateRegular(int sides, double radius)
        {
            if (sides < 3)
                throw new ArgumentException("Polygon must have at least 3 sides");

            var vertices = new Vector2D[sides];
            double angleStep = 2 * Math.PI / sides;

            for (int i = 0; i < sides; i++)
            {
                double angle = i * angleStep - Math.PI / 2; // Start from top
                vertices[i] = new Vector2D(
                    radius * Math.Cos(angle),
                    radius * Math.Sin(angle)
                );
            }

            return new PolygonShape(vertices);
        }

        /// <summary>
        /// Creates a triangle from three vertices.
        /// </summary>
        public static PolygonShape CreateTriangle(Vector2D a, Vector2D b, Vector2D c)
            => new(new[] { a, b, c });

        private void ComputeNormals()
        {
            for (int i = 0; i < _vertices.Length; i++)
            {
                int next = (i + 1) % _vertices.Length;
                var edge = _vertices[next] - _vertices[i];
                // Outward normal (perpendicular, normalized)
                _normals[i] = new Vector2D(edge.Y, -edge.X).Normalized;
            }
        }

        public override AABB2D ComputeAABB(Vector2D position, double rotation)
        {
            if (_vertices.Length == 0)
                return new AABB2D(position, position);

            var first = position + Vector2D.Rotate(_vertices[0], rotation);
            var min = first;
            var max = first;

            for (int i = 1; i < _vertices.Length; i++)
            {
                var v = position + Vector2D.Rotate(_vertices[i], rotation);
                min = Vector2D.Min(min, v);
                max = Vector2D.Max(max, v);
            }

            return new AABB2D(min, max);
        }

        public override MassData2D ComputeMass(double density)
        {
            if (_vertices.Length < 3)
                return new MassData2D();

            // Compute polygon centroid, area, and moment of inertia
            double area = 0;
            Vector2D centroid = Vector2D.Zero;
            double inertia = 0;

            for (int i = 0; i < _vertices.Length; i++)
            {
                var v1 = _vertices[i];
                var v2 = _vertices[(i + 1) % _vertices.Length];

                double cross = Vector2D.Cross(v1, v2);
                double triangleArea = cross * 0.5;
                area += triangleArea;

                // Centroid contribution
                centroid += (v1 + v2) * triangleArea;

                // Inertia contribution (using triangular decomposition)
                double intX2 = v1.X * v1.X + v1.X * v2.X + v2.X * v2.X;
                double intY2 = v1.Y * v1.Y + v1.Y * v2.Y + v2.Y * v2.Y;
                inertia += (0.25 / 3.0) * cross * (intX2 + intY2);
            }

            area = Math.Abs(area);
            double mass = density * area;

            if (area > PhysicsConstants.Epsilon)
            {
                centroid /= (3.0 * area);
            }

            inertia = Math.Abs(inertia) * density;

            // Shift to center of mass using parallel axis theorem
            inertia -= mass * centroid.MagnitudeSquared;

            return new MassData2D
            {
                Mass = mass,
                Inertia = Math.Abs(inertia),
                Center = centroid
            };
        }
    }

    /// <summary>
    /// Capsule shape (line segment with radius) for 2D physics.
    /// </summary>
    public class CapsuleShape : Shape2D
    {
        /// <summary>Half-length of the capsule's central segment.</summary>
        public double HalfLength { get; set; } = 0.5;

        /// <summary>Radius of the capsule.</summary>
        public double Radius { get; set; } = 0.25;

        /// <summary>Whether the capsule is vertical (true) or horizontal (false).</summary>
        public bool IsVertical { get; set; } = true;

        public override ShapeType2D Type => ShapeType2D.Capsule;

        public CapsuleShape() { }

        public CapsuleShape(double halfLength, double radius, bool vertical = true)
        {
            HalfLength = halfLength;
            Radius = radius;
            IsVertical = vertical;
        }

        /// <summary>
        /// Gets the two endpoints of the capsule's central segment in local space.
        /// </summary>
        public (Vector2D A, Vector2D B) GetEndpoints()
        {
            if (IsVertical)
            {
                return (new Vector2D(0, -HalfLength), new Vector2D(0, HalfLength));
            }
            return (new Vector2D(-HalfLength, 0), new Vector2D(HalfLength, 0));
        }

        public override AABB2D ComputeAABB(Vector2D position, double rotation)
        {
            var (a, b) = GetEndpoints();
            a = position + Vector2D.Rotate(a, rotation);
            b = position + Vector2D.Rotate(b, rotation);

            var min = Vector2D.Min(a, b) - new Vector2D(Radius);
            var max = Vector2D.Max(a, b) + new Vector2D(Radius);

            return new AABB2D(min, max);
        }

        public override MassData2D ComputeMass(double density)
        {
            // Capsule = rectangle + 2 half-circles = rectangle + 1 circle
            double rectArea = HalfLength * 2 * Radius * 2;
            double circleArea = Math.PI * Radius * Radius;
            double totalArea = rectArea + circleArea;
            double mass = density * totalArea;

            // Inertia for rectangle around center
            double rectI = density * rectArea * (
                (HalfLength * 2) * (HalfLength * 2) / 12.0 +
                (Radius * 2) * (Radius * 2) / 12.0
            );

            // Inertia for the two half-circles (as one circle, then shifted)
            double circleI = density * circleArea * Radius * Radius * 0.5;
            // Parallel axis theorem for the circle parts
            circleI += density * circleArea * HalfLength * HalfLength;

            return new MassData2D
            {
                Mass = mass,
                Inertia = rectI + circleI,
                Center = Vector2D.Zero
            };
        }
    }

    /// <summary>
    /// Edge (line segment) shape for 2D physics. Typically used for static geometry.
    /// </summary>
    public class EdgeShape : Shape2D
    {
        /// <summary>Start point of the edge.</summary>
        public Vector2D Start { get; set; }

        /// <summary>End point of the edge.</summary>
        public Vector2D End { get; set; }

        /// <summary>Optional ghost vertex for smooth chain collision.</summary>
        public Vector2D? GhostStart { get; set; }

        /// <summary>Optional ghost vertex for smooth chain collision.</summary>
        public Vector2D? GhostEnd { get; set; }

        public override ShapeType2D Type => ShapeType2D.Edge;

        public EdgeShape() { }

        public EdgeShape(Vector2D start, Vector2D end)
        {
            Start = start;
            End = end;
        }

        public Vector2D Direction => (End - Start).Normalized;
        public Vector2D Normal => Direction.Perpendicular;
        public double Length => Vector2D.Distance(Start, End);

        public override AABB2D ComputeAABB(Vector2D position, double rotation)
        {
            var a = position + Vector2D.Rotate(Start, rotation);
            var b = position + Vector2D.Rotate(End, rotation);

            var min = Vector2D.Min(a, b);
            var max = Vector2D.Max(a, b);

            // Add small margin for numerical stability
            var margin = new Vector2D(PhysicsConstants.Epsilon * 10);
            return new AABB2D(min - margin, max + margin);
        }

        public override MassData2D ComputeMass(double density)
        {
            // Edges are typically used for static bodies, so mass is 0
            return new MassData2D
            {
                Mass = 0,
                Inertia = 0,
                Center = (Start + End) * 0.5
            };
        }
    }

    /// <summary>
    /// Chain shape (series of connected edges) for 2D physics.
    /// </summary>
    public class ChainShape : Shape2D
    {
        private Vector2D[] _vertices;

        /// <summary>Gets the vertices of the chain.</summary>
        public Vector2D[] Vertices => _vertices;

        /// <summary>Whether the chain forms a closed loop.</summary>
        public bool IsLoop { get; private set; }

        public int VertexCount => _vertices?.Length ?? 0;
        public int EdgeCount => IsLoop ? VertexCount : Math.Max(0, VertexCount - 1);

        public override ShapeType2D Type => ShapeType2D.Chain;

        public ChainShape()
        {
            _vertices = Array.Empty<Vector2D>();
        }

        public ChainShape(Vector2D[] vertices, bool loop = false)
        {
            _vertices = Array.Empty<Vector2D>();
            if (loop)
                CreateLoop(vertices);
            else
                CreateChain(vertices);
        }

        /// <summary>
        /// Creates an open chain from the given vertices.
        /// </summary>
        public void CreateChain(Vector2D[] vertices)
        {
            if (vertices.Length < 2)
                throw new ArgumentException("Chain must have at least 2 vertices");

            _vertices = new Vector2D[vertices.Length];
            Array.Copy(vertices, _vertices, vertices.Length);
            IsLoop = false;
        }

        /// <summary>
        /// Creates a closed loop from the given vertices.
        /// </summary>
        public void CreateLoop(Vector2D[] vertices)
        {
            if (vertices.Length < 3)
                throw new ArgumentException("Loop must have at least 3 vertices");

            _vertices = new Vector2D[vertices.Length];
            Array.Copy(vertices, _vertices, vertices.Length);
            IsLoop = true;
        }

        /// <summary>
        /// Gets an edge shape for the given edge index.
        /// </summary>
        public EdgeShape GetEdge(int index)
        {
            if (index < 0 || index >= EdgeCount)
                throw new ArgumentOutOfRangeException(nameof(index));

            int next = (index + 1) % VertexCount;
            var edge = new EdgeShape(_vertices[index], _vertices[next]);

            // Set ghost vertices for smooth collision
            if (index > 0 || IsLoop)
            {
                int prev = (index - 1 + VertexCount) % VertexCount;
                edge.GhostStart = _vertices[prev];
            }

            if (index < EdgeCount - 1 || IsLoop)
            {
                int nextNext = (next + 1) % VertexCount;
                edge.GhostEnd = _vertices[nextNext];
            }

            return edge;
        }

        public override AABB2D ComputeAABB(Vector2D position, double rotation)
        {
            if (_vertices.Length == 0)
                return new AABB2D(position, position);

            var first = position + Vector2D.Rotate(_vertices[0], rotation);
            var min = first;
            var max = first;

            for (int i = 1; i < _vertices.Length; i++)
            {
                var v = position + Vector2D.Rotate(_vertices[i], rotation);
                min = Vector2D.Min(min, v);
                max = Vector2D.Max(max, v);
            }

            return new AABB2D(min, max);
        }

        public override MassData2D ComputeMass(double density)
        {
            // Chains are typically used for static bodies
            return new MassData2D
            {
                Mass = 0,
                Inertia = 0,
                Center = Vector2D.Zero
            };
        }
    }
}
