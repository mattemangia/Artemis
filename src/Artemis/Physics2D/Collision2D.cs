using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Contact point between two colliding bodies.
    /// </summary>
    public struct Contact2D
    {
        /// <summary>Contact point in world space.</summary>
        public Vector2D Point;

        /// <summary>Contact normal (points from body A to body B).</summary>
        public Vector2D Normal;

        /// <summary>Penetration depth (positive = overlapping).</summary>
        public double Depth;

        /// <summary>Accumulated normal impulse for warm starting.</summary>
        public double NormalImpulse;

        /// <summary>Accumulated tangent impulse for warm starting.</summary>
        public double TangentImpulse;

        // Solver state, valid for the duration of a step.
        internal Vector2D LocalA;
        internal Vector2D LocalB;
        internal Vector2D RA;
        internal Vector2D RB;
        internal double NormalMass;
        internal double TangentMass;
        internal double TargetNormalVelocity;
    }

    /// <summary>
    /// Collision manifold between two bodies.
    /// </summary>
    public class Manifold2D
    {
        /// <summary>First body in collision.</summary>
        public RigidBody2D BodyA { get; set; } = null!;

        /// <summary>Second body in collision.</summary>
        public RigidBody2D BodyB { get; set; } = null!;

        /// <summary>Contact points.</summary>
        public Contact2D[] Contacts { get; } = new Contact2D[2];

        /// <summary>Contact normal (from A to B).</summary>
        public Vector2D Normal { get; set; }

        /// <summary>Penetration depth.</summary>
        public double Penetration { get; set; }

        /// <summary>Number of valid contacts.</summary>
        public int ContactCount { get; set; }

        /// <summary>Combined friction.</summary>
        public double Friction { get; set; }

        /// <summary>Combined restitution.</summary>
        public double Restitution { get; set; }

        /// <summary>Whether this manifold represents a trigger overlap.</summary>
        public bool IsTrigger { get; set; }

        /// <summary>Whether the collision is still active.</summary>
        public bool IsActive { get; set; } = true;
    }

    /// <summary>
    /// 2D collision detection algorithms.
    /// </summary>
    public static class CollisionDetector2D
    {
        #region Circle vs Circle

        /// <summary>
        /// Detects collision between two circles.
        /// </summary>
        public static bool CircleVsCircle(
            Vector2D posA, double radiusA,
            Vector2D posB, double radiusB,
            out Vector2D normal, out double depth)
        {
            normal = Vector2D.Zero;
            depth = 0;

            var d = posB - posA;
            double distSq = d.MagnitudeSquared;
            double radiusSum = radiusA + radiusB;

            if (distSq >= radiusSum * radiusSum)
                return false;

            double dist = Math.Sqrt(distSq);
            if (dist > PhysicsConstants.Epsilon)
            {
                normal = d / dist;
                depth = radiusSum - dist;
            }
            else
            {
                // Circles are at the same position
                normal = Vector2D.Up;
                depth = radiusSum;
            }

            return true;
        }

        /// <summary>
        /// Creates a collision manifold for two circles.
        /// </summary>
        public static bool CircleVsCircle(RigidBody2D a, RigidBody2D b, Manifold2D manifold)
        {
            if (a.Shape is not CircleShape circleA || b.Shape is not CircleShape circleB)
                return false;

            var posA = a.Position + Vector2D.Rotate(circleA.Offset, a.Rotation);
            var posB = b.Position + Vector2D.Rotate(circleB.Offset, b.Rotation);

            if (!CircleVsCircle(posA, circleA.Radius, posB, circleB.Radius,
                out var normal, out var depth))
                return false;

            manifold.BodyA = a;
            manifold.BodyB = b;
            manifold.ContactCount = 1;
            manifold.Contacts[0] = new Contact2D
            {
                Point = posA + normal * circleA.Radius,
                Normal = normal,
                Depth = depth
            };

            return true;
        }

        #endregion

        #region Circle vs Box

        /// <summary>
        /// Detects collision between a circle and a box.
        /// </summary>
        public static bool CircleVsBox(RigidBody2D circle, RigidBody2D box, Manifold2D manifold)
        {
            if (circle.Shape is not CircleShape circleShape || box.Shape is not BoxShape boxShape)
                return false;

            var circlePos = circle.Position + Vector2D.Rotate(circleShape.Offset, circle.Rotation);

            // Transform circle center to box's local space
            var localCircle = box.WorldToLocal(circlePos) - boxShape.Offset;

            // Find closest point on box to circle
            double closestX = Math.Clamp(localCircle.X, -boxShape.HalfWidth, boxShape.HalfWidth);
            double closestY = Math.Clamp(localCircle.Y, -boxShape.HalfHeight, boxShape.HalfHeight);
            var closest = new Vector2D(closestX, closestY);

            // Check if circle is inside the box
            bool inside = Math.Abs(localCircle.X) <= boxShape.HalfWidth &&
                          Math.Abs(localCircle.Y) <= boxShape.HalfHeight;

            // The manifold normal points from A (the circle) to B (the box), i.e.
            // opposite to the direction that pushes the circle out of the box.
            Vector2D normal;
            double depth;

            if (inside)
            {
                // Centre inside the box: leave through the nearest face.
                double dx = boxShape.HalfWidth - Math.Abs(localCircle.X);
                double dy = boxShape.HalfHeight - Math.Abs(localCircle.Y);

                Vector2D outward;
                if (dx < dy)
                {
                    double sign = localCircle.X >= 0 ? 1 : -1;
                    closest = new Vector2D(sign * boxShape.HalfWidth, localCircle.Y);
                    outward = new Vector2D(sign, 0);
                    depth = circleShape.Radius + dx;
                }
                else
                {
                    double sign = localCircle.Y >= 0 ? 1 : -1;
                    closest = new Vector2D(localCircle.X, sign * boxShape.HalfHeight);
                    outward = new Vector2D(0, sign);
                    depth = circleShape.Radius + dy;
                }

                normal = -Vector2D.Rotate(outward, box.Rotation);
            }
            else
            {
                var d = localCircle - closest;
                double distSq = d.MagnitudeSquared;

                if (distSq >= circleShape.Radius * circleShape.Radius)
                    return false;

                double dist = Math.Sqrt(distSq);
                if (dist < PhysicsConstants.Epsilon)
                    return false;

                normal = -Vector2D.Rotate(d / dist, box.Rotation);
                depth = circleShape.Radius - dist;
            }

            // Transform contact point to world space
            var contactPoint = box.LocalToWorld(closest + boxShape.Offset);

            manifold.BodyA = circle;
            manifold.BodyB = box;
            manifold.ContactCount = 1;
            manifold.Contacts[0] = new Contact2D
            {
                Point = contactPoint,
                Normal = normal,
                Depth = depth
            };

            return true;
        }

        #endregion

        #region Box vs Box (SAT)

        /// <summary>
        /// Detects collision between two boxes using Separating Axis Theorem.
        /// </summary>
        public static bool BoxVsBox(RigidBody2D a, RigidBody2D b, Manifold2D manifold)
        {
            if (a.Shape is not BoxShape boxA || b.Shape is not BoxShape boxB)
                return false;

            var verticesA = boxA.GetWorldVertices(a.Position, a.Rotation);
            var verticesB = boxB.GetWorldVertices(b.Position, b.Rotation);

            // Get axes to test (4 axes: 2 from each box)
            var axes = new Vector2D[4];
            for (int i = 0; i < 2; i++)
            {
                var edge = verticesA[(i + 1) % 4] - verticesA[i];
                axes[i] = new Vector2D(-edge.Y, edge.X).Normalized;
            }
            for (int i = 0; i < 2; i++)
            {
                var edge = verticesB[(i + 1) % 4] - verticesB[i];
                axes[i + 2] = new Vector2D(-edge.Y, edge.X).Normalized;
            }

            double minOverlap = double.MaxValue;
            Vector2D minAxis = Vector2D.Zero;
            int bestAxisIndex = 0;
            int currentAxisIndex = 0;

            foreach (var axis in axes)
            {
                var (minA, maxA) = ProjectPolygon(verticesA, axis);
                var (minB, maxB) = ProjectPolygon(verticesB, axis);

                double overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);
                if (overlap <= 0)
                    return false;

                if (overlap < minOverlap)
                {
                    minOverlap = overlap;
                    minAxis = axis;
                    bestAxisIndex = currentAxisIndex;
                }
                currentAxisIndex++;
            }

            // Ensure normal points from A to B
            var centerDiff = b.Position - a.Position;
            if (Vector2D.Dot(centerDiff, minAxis) < 0)
                minAxis = -minAxis;

            // Find contact points
            FindBoxContactPoints(verticesA, verticesB, minAxis, manifold, bestAxisIndex);

            manifold.BodyA = a;
            manifold.BodyB = b;

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                var contact = manifold.Contacts[i];
                contact.Normal = minAxis;
                contact.Depth = minOverlap;
                manifold.Contacts[i] = contact;
            }

            return true;
        }

        private static (double min, double max) ProjectPolygon(Vector2D[] vertices, Vector2D axis)
        {
            double min = Vector2D.Dot(vertices[0], axis);
            double max = min;

            for (int i = 1; i < vertices.Length; i++)
            {
                double proj = Vector2D.Dot(vertices[i], axis);
                if (proj < min) min = proj;
                if (proj > max) max = proj;
            }

            return (min, max);
        }

        private static void FindBoxContactPoints(Vector2D[] verticesA, Vector2D[] verticesB,
            Vector2D normal, Manifold2D manifold, int bestAxisIndex)
        {
            // Axes 0-1 are from A, 2-3 are from B
            if (bestAxisIndex < 2)
            {
                // Reference face on A, Incident face on B
                // Find support point on B (deepest in A's direction, i.e., against normal)
                // Note: Normal points A -> B. So we want point on B furthest in -Normal direction.
                int supportB = FindSupportPoint(verticesB, -normal);
                manifold.ContactCount = 1;
                manifold.Contacts[0] = new Contact2D
                {
                    Point = verticesB[supportB]
                };
            }
            else
            {
                // Reference face on B, Incident face on A
                // Find support point on A (deepest in B's direction, i.e., along normal)
                int supportA = FindSupportPoint(verticesA, normal);
                manifold.ContactCount = 1;
                manifold.Contacts[0] = new Contact2D
                {
                    Point = verticesA[supportA]
                };
            }
        }

        private static int FindSupportPoint(Vector2D[] vertices, Vector2D direction)
        {
            int bestIndex = 0;
            double bestProj = Vector2D.Dot(vertices[0], direction);

            for (int i = 1; i < vertices.Length; i++)
            {
                double proj = Vector2D.Dot(vertices[i], direction);
                if (proj > bestProj)
                {
                    bestProj = proj;
                    bestIndex = i;
                }
            }

            return bestIndex;
        }

        #endregion

        #region Polygon vs Polygon (SAT)

        /// <summary>
        /// Detects collision between two convex polygons using SAT.
        /// </summary>
        public static bool PolygonVsPolygon(RigidBody2D a, RigidBody2D b, Manifold2D manifold)
        {
            if (a.Shape is not PolygonShape polyA || b.Shape is not PolygonShape polyB)
                return false;

            // Transform vertices to world space
            var verticesA = new Vector2D[polyA.VertexCount];
            var verticesB = new Vector2D[polyB.VertexCount];

            for (int i = 0; i < polyA.VertexCount; i++)
                verticesA[i] = a.LocalToWorld(polyA.Vertices[i]);

            for (int i = 0; i < polyB.VertexCount; i++)
                verticesB[i] = b.LocalToWorld(polyB.Vertices[i]);

            // Test all axes from both polygons
            double minOverlap = double.MaxValue;
            Vector2D minAxis = Vector2D.Zero;
            bool minAxisFromA = true;

            // Test axes from polygon A
            for (int i = 0; i < polyA.VertexCount; i++)
            {
                var axis = Vector2D.Rotate(polyA.Normals[i], a.Rotation);
                var (minA, maxA) = ProjectPolygon(verticesA, axis);
                var (minB, maxB) = ProjectPolygon(verticesB, axis);

                double overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);
                if (overlap <= 0)
                    return false;

                if (overlap < minOverlap)
                {
                    minOverlap = overlap;
                    minAxis = axis;
                    minAxisFromA = true;
                }
            }

            // Test axes from polygon B
            for (int i = 0; i < polyB.VertexCount; i++)
            {
                var axis = Vector2D.Rotate(polyB.Normals[i], b.Rotation);
                var (minA, maxA) = ProjectPolygon(verticesA, axis);
                var (minB, maxB) = ProjectPolygon(verticesB, axis);

                double overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);
                if (overlap <= 0)
                    return false;

                if (overlap < minOverlap)
                {
                    minOverlap = overlap;
                    minAxis = axis;
                    minAxisFromA = false;
                }
            }

            // Ensure normal points from A to B
            var centerDiff = b.Position - a.Position;
            if (Vector2D.Dot(centerDiff, minAxis) < 0)
                minAxis = -minAxis;

            manifold.BodyA = a;
            manifold.BodyB = b;
            manifold.ContactCount = 1;

            if (minAxisFromA)
            {
                // A provides normal. B is Incident.
                int supportB = FindSupportPoint(verticesB, -minAxis);
                manifold.Contacts[0] = new Contact2D
                {
                    Point = verticesB[supportB],
                    Normal = minAxis,
                    Depth = minOverlap
                };
            }
            else
            {
                // B provides normal. A is Incident.
                int supportA = FindSupportPoint(verticesA, minAxis);
                manifold.Contacts[0] = new Contact2D
                {
                    Point = verticesA[supportA],
                    Normal = minAxis,
                    Depth = minOverlap
                };
            }

            return true;
        }

        #endregion

        #region Circle vs Edge

        /// <summary>
        /// Detects collision between a circle and an edge.
        /// </summary>
        public static bool CircleVsEdge(RigidBody2D circle, RigidBody2D edge, Manifold2D manifold)
        {
            if (circle.Shape is not CircleShape circleShape || edge.Shape is not EdgeShape edgeShape)
                return false;

            var circlePos = circle.Position + Vector2D.Rotate(circleShape.Offset, circle.Rotation);
            var start = edge.LocalToWorld(edgeShape.Start);
            var end = edge.LocalToWorld(edgeShape.End);

            // Find closest point on edge to circle center
            var edgeVec = end - start;
            var toCircle = circlePos - start;
            double t = Math.Clamp(Vector2D.Dot(toCircle, edgeVec) / edgeVec.MagnitudeSquared, 0, 1);
            var closest = start + edgeVec * t;

            var d = circlePos - closest;
            double distSq = d.MagnitudeSquared;

            if (distSq >= circleShape.Radius * circleShape.Radius)
                return false;

            double dist = Math.Sqrt(distSq);
            Vector2D normal = dist > PhysicsConstants.Epsilon ? d / dist : edgeShape.Normal;

            manifold.BodyA = circle;
            manifold.BodyB = edge;
            manifold.ContactCount = 1;
            manifold.Contacts[0] = new Contact2D
            {
                Point = closest,
                Normal = normal,
                Depth = circleShape.Radius - dist
            };

            return true;
        }

        #endregion

        #region Point Tests

        /// <summary>
        /// Tests if a point is inside a circle.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool PointInCircle(Vector2D point, Vector2D center, double radius)
            => Vector2D.DistanceSquared(point, center) <= radius * radius;

        /// <summary>
        /// Tests if a point is inside an AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool PointInAABB(Vector2D point, AABB2D aabb)
            => aabb.Contains(point);

        /// <summary>
        /// Tests if a point is inside a convex polygon.
        /// </summary>
        public static bool PointInPolygon(Vector2D point, Vector2D[] vertices)
        {
            int count = vertices.Length;
            bool inside = false;

            for (int i = 0, j = count - 1; i < count; j = i++)
            {
                if ((vertices[i].Y > point.Y) != (vertices[j].Y > point.Y) &&
                    point.X < (vertices[j].X - vertices[i].X) * (point.Y - vertices[i].Y) /
                             (vertices[j].Y - vertices[i].Y) + vertices[i].X)
                {
                    inside = !inside;
                }
            }

            return inside;
        }

        #endregion

        #region Ray Casting

        /// <summary>
        /// Result of a 2D raycast.
        /// </summary>
        public struct RaycastHit2D
        {
            /// <summary>Whether the ray hit something.</summary>
            public bool Hit;

            /// <summary>The body that was hit.</summary>
            public RigidBody2D? Body;

            /// <summary>Hit point in world space.</summary>
            public Vector2D Point;

            /// <summary>Surface normal at hit point.</summary>
            public Vector2D Normal;

            /// <summary>Distance along ray to hit point (0-1 for normalized).</summary>
            public double Fraction;
        }

        /// <summary>
        /// Casts a ray against a circle.
        /// </summary>
        public static bool RaycastCircle(Vector2D origin, Vector2D direction, double maxDistance,
            Vector2D center, double radius, out RaycastHit2D hit)
        {
            hit = default;

            var toCenter = center - origin;
            double b = Vector2D.Dot(toCenter, direction);
            double c = toCenter.MagnitudeSquared - radius * radius;

            // Ray starts inside or misses
            if (c > 0 && b < 0)
                return false;

            double discriminant = b * b - c;
            if (discriminant < 0)
                return false;

            double t = b - Math.Sqrt(discriminant);
            if (t < 0) t = 0; // Inside circle
            if (t > maxDistance)
                return false;

            hit.Hit = true;
            hit.Fraction = t / maxDistance;
            hit.Point = origin + direction * t;
            hit.Normal = (hit.Point - center).Normalized;

            return true;
        }

        /// <summary>
        /// Casts a ray against an AABB.
        /// </summary>
        public static bool RaycastAABB(Vector2D origin, Vector2D direction, double maxDistance,
            AABB2D aabb, out RaycastHit2D hit)
        {
            hit = default;

            double tMin = 0;
            double tMax = maxDistance;

            // X slab
            if (Math.Abs(direction.X) < PhysicsConstants.Epsilon)
            {
                if (origin.X < aabb.Min.X || origin.X > aabb.Max.X)
                    return false;
            }
            else
            {
                double invD = 1.0 / direction.X;
                double t1 = (aabb.Min.X - origin.X) * invD;
                double t2 = (aabb.Max.X - origin.X) * invD;

                if (t1 > t2) (t1, t2) = (t2, t1);

                tMin = Math.Max(tMin, t1);
                tMax = Math.Min(tMax, t2);

                if (tMin > tMax)
                    return false;
            }

            // Y slab
            if (Math.Abs(direction.Y) < PhysicsConstants.Epsilon)
            {
                if (origin.Y < aabb.Min.Y || origin.Y > aabb.Max.Y)
                    return false;
            }
            else
            {
                double invD = 1.0 / direction.Y;
                double t1 = (aabb.Min.Y - origin.Y) * invD;
                double t2 = (aabb.Max.Y - origin.Y) * invD;

                if (t1 > t2) (t1, t2) = (t2, t1);

                tMin = Math.Max(tMin, t1);
                tMax = Math.Min(tMax, t2);

                if (tMin > tMax)
                    return false;
            }

            hit.Hit = true;
            hit.Fraction = tMin / maxDistance;
            hit.Point = origin + direction * tMin;

            // Determine normal
            var localPoint = hit.Point - aabb.Center;
            var absPoint = new Vector2D(Math.Abs(localPoint.X), Math.Abs(localPoint.Y));

            if (absPoint.X > absPoint.Y)
                hit.Normal = new Vector2D(Math.Sign(localPoint.X), 0);
            else
                hit.Normal = new Vector2D(0, Math.Sign(localPoint.Y));

            return true;
        }

        #endregion

        #region Convex Polygons (SAT + clipping)

        /// <summary>
        /// World-space vertices of a convex shape (box, polygon or edge), in
        /// counter-clockwise order. Returns <c>null</c> for other shapes.
        /// </summary>
        public static Vector2D[]? GetConvexWorldVertices(RigidBody2D body)
        {
            Vector2D[] vertices;
            switch (body.Shape)
            {
                case BoxShape box:
                    vertices = box.GetWorldVertices(body.Position, body.Rotation);
                    break;

                case PolygonShape polygon when polygon.VertexCount >= 3:
                    vertices = new Vector2D[polygon.VertexCount];
                    for (int i = 0; i < vertices.Length; i++)
                        vertices[i] = body.LocalToWorld(polygon.Vertices[i]);
                    break;

                case EdgeShape edge:
                    return new[] { body.LocalToWorld(edge.Start), body.LocalToWorld(edge.End) };

                default:
                    return null;
            }

            // Outward normals below assume counter-clockwise winding.
            double area = 0;
            for (int i = 0; i < vertices.Length; i++)
                area += Vector2D.Cross(vertices[i], vertices[(i + 1) % vertices.Length]);
            if (area < 0)
                Array.Reverse(vertices);

            return vertices;
        }

        private static Vector2D EdgeNormal(Vector2D[] vertices, int index)
        {
            var edge = vertices[(index + 1) % vertices.Length] - vertices[index];
            return new Vector2D(edge.Y, -edge.X).Normalized;
        }

        /// <summary>
        /// Largest separation of <paramref name="other"/> along the face normals of
        /// <paramref name="reference"/>. Positive means the shapes do not overlap.
        /// </summary>
        private static double FindMaxSeparation(Vector2D[] reference, Vector2D[] other, out int bestEdge)
        {
            double best = double.NegativeInfinity;
            bestEdge = 0;

            for (int i = 0; i < reference.Length; i++)
            {
                var normal = EdgeNormal(reference, i);
                double separation = double.PositiveInfinity;
                for (int j = 0; j < other.Length; j++)
                    separation = Math.Min(separation, Vector2D.Dot(normal, other[j] - reference[i]));

                if (separation > best)
                {
                    best = separation;
                    bestEdge = i;
                }
            }

            return best;
        }

        /// <summary>
        /// Collision between two convex shapes given by their world-space vertices.
        /// Produces up to two contact points by clipping the incident edge against the
        /// reference face, which is what keeps boxes from rocking on a single corner.
        /// </summary>
        public static bool ConvexVsConvex(RigidBody2D a, Vector2D[] verticesA, RigidBody2D b, Vector2D[] verticesB,
            Manifold2D manifold)
        {
            double separationA = FindMaxSeparation(verticesA, verticesB, out int edgeA);
            if (separationA > 0)
                return false;

            double separationB = FindMaxSeparation(verticesB, verticesA, out int edgeB);
            if (separationB > 0)
                return false;

            // Prefer A as reference unless B's face is clearly the shallower one, so the
            // choice does not flip back and forth between steps.
            bool referenceIsA = separationA >= separationB - 0.0005;

            var reference = referenceIsA ? verticesA : verticesB;
            var incident = referenceIsA ? verticesB : verticesA;
            int refEdge = referenceIsA ? edgeA : edgeB;

            var refNormal = EdgeNormal(reference, refEdge);
            var ref1 = reference[refEdge];
            var ref2 = reference[(refEdge + 1) % reference.Length];

            // Incident edge: the one whose normal opposes the reference normal most.
            int incEdge = 0;
            double minDot = double.PositiveInfinity;
            for (int i = 0; i < incident.Length; i++)
            {
                double dot = Vector2D.Dot(EdgeNormal(incident, i), refNormal);
                if (dot < minDot)
                {
                    minDot = dot;
                    incEdge = i;
                }
            }

            var inc1 = incident[incEdge];
            var inc2 = incident[(incEdge + 1) % incident.Length];

            // Clip the incident edge to the side planes of the reference edge.
            var tangent = (ref2 - ref1).Normalized;
            if (!ClipSegment(ref inc1, ref inc2, -tangent, -Vector2D.Dot(tangent, ref1)) ||
                !ClipSegment(ref inc1, ref inc2, tangent, Vector2D.Dot(tangent, ref2)))
                return false;

            var normal = referenceIsA ? refNormal : -refNormal;

            int count = 0;
            double deepest = 0;
            AddClippedPoint(manifold, ref count, ref deepest, inc1, refNormal, ref1, normal);
            if ((inc2 - inc1).MagnitudeSquared > 1e-12)
                AddClippedPoint(manifold, ref count, ref deepest, inc2, refNormal, ref1, normal);

            if (count == 0)
                return false;

            manifold.BodyA = a;
            manifold.BodyB = b;
            manifold.ContactCount = count;
            manifold.Normal = normal;
            manifold.Penetration = deepest;
            return true;
        }

        private static void AddClippedPoint(Manifold2D manifold, ref int count, ref double deepest,
            Vector2D point, Vector2D refNormal, Vector2D refPoint, Vector2D normal)
        {
            double separation = Vector2D.Dot(refNormal, point - refPoint);
            if (separation > 0)
                return;

            manifold.Contacts[count++] = new Contact2D { Point = point, Normal = normal, Depth = -separation };
            deepest = Math.Max(deepest, -separation);
        }

        /// <summary>
        /// Clips a segment to the half plane <c>dot(normal, p) &lt;= offset</c>.
        /// Returns false if nothing is left.
        /// </summary>
        private static bool ClipSegment(ref Vector2D p1, ref Vector2D p2, Vector2D normal, double offset)
        {
            double d1 = Vector2D.Dot(normal, p1) - offset;
            double d2 = Vector2D.Dot(normal, p2) - offset;

            if (d1 > 0 && d2 > 0)
                return false;

            if (d1 > 0)
                p1 += (p2 - p1) * (d1 / (d1 - d2));
            else if (d2 > 0)
                p2 += (p1 - p2) * (d2 / (d2 - d1));

            return true;
        }

        /// <summary>
        /// Collision between a circle (body A) and a convex polygon or edge (body B).
        /// </summary>
        public static bool CircleVsConvex(RigidBody2D circle, RigidBody2D convex, Vector2D[] vertices, Manifold2D manifold)
        {
            if (circle.Shape is not CircleShape circleShape)
                return false;

            var center = circle.Position + Vector2D.Rotate(circleShape.Offset, circle.Rotation);
            double radius = circleShape.Radius;

            // Closest point on the boundary, and whether the centre is inside.
            double bestDistanceSq = double.PositiveInfinity;
            var closest = vertices[0];
            double maxSeparation = double.NegativeInfinity;
            int deepestEdge = 0;
            int edgeCount = vertices.Length == 2 ? 1 : vertices.Length;

            for (int i = 0; i < edgeCount; i++)
            {
                var v1 = vertices[i];
                var v2 = vertices[(i + 1) % vertices.Length];
                var edge = v2 - v1;
                double lengthSq = edge.MagnitudeSquared;
                double t = lengthSq > PhysicsConstants.Epsilon
                    ? Math.Clamp(Vector2D.Dot(center - v1, edge) / lengthSq, 0, 1)
                    : 0;
                var point = v1 + edge * t;
                double distanceSq = (center - point).MagnitudeSquared;
                if (distanceSq < bestDistanceSq)
                {
                    bestDistanceSq = distanceSq;
                    closest = point;
                }

                double separation = Vector2D.Dot(EdgeNormal(vertices, i), center - v1);
                if (separation > maxSeparation)
                {
                    maxSeparation = separation;
                    deepestEdge = i;
                }
            }

            bool inside = vertices.Length >= 3 && maxSeparation <= 0;

            Vector2D outward; // direction that pushes the circle out of the polygon
            double depth;
            Vector2D contact;

            if (inside)
            {
                outward = EdgeNormal(vertices, deepestEdge);
                depth = radius - maxSeparation;
                contact = center - outward * maxSeparation;
            }
            else
            {
                if (bestDistanceSq >= radius * radius)
                    return false;

                double distance = Math.Sqrt(bestDistanceSq);
                if (distance < PhysicsConstants.Epsilon)
                    outward = EdgeNormal(vertices, deepestEdge);
                else
                    outward = (center - closest) / distance;

                depth = radius - distance;
                contact = closest;
            }

            manifold.BodyA = circle;
            manifold.BodyB = convex;
            manifold.ContactCount = 1;
            manifold.Normal = -outward;
            manifold.Penetration = depth;
            manifold.Contacts[0] = new Contact2D { Point = contact, Normal = -outward, Depth = depth };
            return true;
        }

        /// <summary>
        /// Collides a body against each segment of a chain and keeps the deepest result.
        /// </summary>
        private static bool AnyVsChain(RigidBody2D other, RigidBody2D chainBody, ChainShape chain, Manifold2D manifold)
        {
            var otherVertices = other.Shape is CircleShape ? null : GetConvexWorldVertices(other);
            if (otherVertices == null && other.Shape is not CircleShape)
                return false;

            var scratch = new Manifold2D();
            var segment = new Vector2D[2];
            bool found = false;
            double deepest = 0;

            for (int i = 0; i < chain.EdgeCount; i++)
            {
                segment[0] = chainBody.LocalToWorld(chain.Vertices[i]);
                segment[1] = chainBody.LocalToWorld(chain.Vertices[(i + 1) % chain.VertexCount]);

                bool hit = otherVertices == null
                    ? CircleVsConvex(other, chainBody, segment, scratch)
                    : ConvexVsConvex(other, otherVertices, chainBody, segment, scratch);

                if (!hit || (found && scratch.Penetration <= deepest))
                    continue;

                found = true;
                deepest = scratch.Penetration;
                manifold.BodyA = scratch.BodyA;
                manifold.BodyB = scratch.BodyB;
                manifold.ContactCount = scratch.ContactCount;
                manifold.Normal = scratch.Normal;
                manifold.Penetration = scratch.Penetration;
                for (int c = 0; c < scratch.ContactCount; c++)
                    manifold.Contacts[c] = scratch.Contacts[c];
            }

            return found;
        }

        #endregion

        #region Generic Dispatch

        /// <summary>
        /// Detects collision between two bodies, dispatching to the appropriate algorithm.
        /// </summary>
        public static bool DetectCollision(RigidBody2D a, RigidBody2D b, Manifold2D manifold)
        {
            if (a.Shape == null || b.Shape == null)
                return false;

            // AABB early-out
            if (!a.AABB.Overlaps(b.AABB))
                return false;

            // Dispatch based on shape types
            var typeA = a.Shape.Type;
            var typeB = b.Shape.Type;

            // Ensure consistent ordering for asymmetric pairs
            if (typeA > typeB)
            {
                (a, b) = (b, a);
                (typeA, typeB) = (typeB, typeA);
            }

            bool result;
            if (typeB == ShapeType2D.Chain)
            {
                result = typeA != ShapeType2D.Chain && AnyVsChain(a, b, (ChainShape)b.Shape!, manifold);
            }
            else if (typeA == ShapeType2D.Circle)
            {
                if (typeB == ShapeType2D.Circle)
                    result = CircleVsCircle(a, b, manifold);
                else if (typeB == ShapeType2D.Box)
                    result = CircleVsBox(a, b, manifold);
                else
                    result = GetConvexWorldVertices(b) is { } vertices && CircleVsConvex(a, b, vertices, manifold);
            }
            else
            {
                // Boxes, polygons and edges in any combination.
                var verticesA = GetConvexWorldVertices(a);
                var verticesB = GetConvexWorldVertices(b);
                result = verticesA != null && verticesB != null &&
                         !(verticesA.Length == 2 && verticesB.Length == 2) &&
                         ConvexVsConvex(a, verticesA, b, verticesB, manifold);
            }

            if (result)
            {
                // Calculate combined material properties
                manifold.Friction = Math.Sqrt(a.Material.StaticFriction * b.Material.StaticFriction);
                manifold.Restitution = Math.Max(a.Material.Restitution, b.Material.Restitution);
            }

            return result;
        }

        #endregion
    }

    /// <summary>
    /// Resolves collisions by applying impulses.
    /// </summary>
    public static class CollisionResolver2D
    {
        /// <summary>
        /// Resolves a collision manifold by applying impulses.
        /// </summary>
        public static void ResolveCollision(Manifold2D manifold, int iterations = 1)
        {
            // One-shot resolution of a single manifold. The physics world instead keeps
            // manifolds alive across steps and warm starts them (see ContactSolver2D).
            ContactSolver2D.PreStep(manifold, restitutionThreshold: 0, warmStart: false);

            for (int iter = 0; iter < Math.Max(iterations, 1); iter++)
                ContactSolver2D.SolveVelocity(manifold);
        }

        /// <summary>
        /// Corrects penetration by directly moving bodies apart.
        /// </summary>
        public static void CorrectPositions(Manifold2D manifold, double slop = 0.01, double percent = 0.4)
        {
            var a = manifold.BodyA;
            var b = manifold.BodyB;

            if (manifold.ContactCount == 0)
                return;

            double maxDepth = 0;
            Vector2D normal = Vector2D.Zero;

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                if (manifold.Contacts[i].Depth > maxDepth)
                {
                    maxDepth = manifold.Contacts[i].Depth;
                    normal = manifold.Contacts[i].Normal;
                }
            }

            double correction = Math.Max(maxDepth - slop, 0) / (a.InverseMass + b.InverseMass) * percent;

            if (correction > PhysicsConstants.Epsilon)
            {
                var correctionVec = normal * correction;

                if (a.BodyType == BodyType2D.Dynamic)
                    a.Position -= correctionVec * a.InverseMass;

                if (b.BodyType == BodyType2D.Dynamic)
                    b.Position += correctionVec * b.InverseMass;
            }
        }
    }
}
