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

        /// <summary>Number of valid contacts.</summary>
        public int ContactCount { get; set; }

        /// <summary>Combined friction.</summary>
        public double Friction { get; set; }

        /// <summary>Combined restitution.</summary>
        public double Restitution { get; set; }

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
            var localCircle = box.WorldToLocal(circlePos);

            // Find closest point on box to circle
            double closestX = Math.Clamp(localCircle.X, -boxShape.HalfWidth, boxShape.HalfWidth);
            double closestY = Math.Clamp(localCircle.Y, -boxShape.HalfHeight, boxShape.HalfHeight);
            var closest = new Vector2D(closestX, closestY);

            // Check if circle is inside the box
            bool inside = Math.Abs(localCircle.X) < boxShape.HalfWidth &&
                          Math.Abs(localCircle.Y) < boxShape.HalfHeight;

            Vector2D normal;
            double depth;

            if (inside)
            {
                // Find which edge is closest
                double dx = boxShape.HalfWidth - Math.Abs(localCircle.X);
                double dy = boxShape.HalfHeight - Math.Abs(localCircle.Y);

                if (dx < dy)
                {
                    closest = new Vector2D(
                        localCircle.X > 0 ? boxShape.HalfWidth : -boxShape.HalfWidth,
                        localCircle.Y
                    );
                }
                else
                {
                    closest = new Vector2D(
                        localCircle.X,
                        localCircle.Y > 0 ? boxShape.HalfHeight : -boxShape.HalfHeight
                    );
                }

                var d = localCircle - closest;
                double dist = d.Magnitude;
                normal = Vector2D.Rotate(d.Normalized, box.Rotation);
                depth = circleShape.Radius + dist;
            }
            else
            {
                var d = localCircle - closest;
                double distSq = d.MagnitudeSquared;

                if (distSq >= circleShape.Radius * circleShape.Radius)
                    return false;

                double dist = Math.Sqrt(distSq);
                normal = Vector2D.Rotate(d / dist, box.Rotation);
                depth = circleShape.Radius - dist;
            }

            // Transform contact point to world space
            var contactPoint = box.LocalToWorld(closest);

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
                }
            }

            // Ensure normal points from A to B
            var centerDiff = b.Position - a.Position;
            if (Vector2D.Dot(centerDiff, minAxis) < 0)
                minAxis = -minAxis;

            // Find contact points
            FindBoxContactPoints(verticesA, verticesB, minAxis, manifold);

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
            Vector2D normal, Manifold2D manifold)
        {
            // Find the most penetrating vertex from each polygon
            int supportA = FindSupportPoint(verticesA, -normal);
            int supportB = FindSupportPoint(verticesB, normal);

            // Simple contact: use support points
            manifold.ContactCount = 1;
            manifold.Contacts[0] = new Contact2D
            {
                Point = (verticesA[supportA] + verticesB[supportB]) * 0.5
            };
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
                }
            }

            // Ensure normal points from A to B
            var centerDiff = b.Position - a.Position;
            if (Vector2D.Dot(centerDiff, minAxis) < 0)
                minAxis = -minAxis;

            // Find contact point (simplified: use support point)
            int supportB = FindSupportPoint(verticesB, -minAxis);

            manifold.BodyA = a;
            manifold.BodyB = b;
            manifold.ContactCount = 1;
            manifold.Contacts[0] = new Contact2D
            {
                Point = verticesB[supportB],
                Normal = minAxis,
                Depth = minOverlap
            };

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

            bool result = (typeA, typeB) switch
            {
                (ShapeType2D.Circle, ShapeType2D.Circle) => CircleVsCircle(a, b, manifold),
                (ShapeType2D.Circle, ShapeType2D.Box) => CircleVsBox(a, b, manifold),
                (ShapeType2D.Circle, ShapeType2D.Edge) => CircleVsEdge(a, b, manifold),
                (ShapeType2D.Box, ShapeType2D.Box) => BoxVsBox(a, b, manifold),
                (ShapeType2D.Polygon, ShapeType2D.Polygon) => PolygonVsPolygon(a, b, manifold),
                _ => false
            };

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
            var a = manifold.BodyA;
            var b = manifold.BodyB;

            for (int iter = 0; iter < iterations; iter++)
            {
                for (int i = 0; i < manifold.ContactCount; i++)
                {
                    ref var contact = ref manifold.Contacts[i];

                    var ra = contact.Point - a.Position;
                    var rb = contact.Point - b.Position;

                    // Relative velocity at contact
                    var va = a.Velocity + Vector2D.Cross(a.AngularVelocity, ra);
                    var vb = b.Velocity + Vector2D.Cross(b.AngularVelocity, rb);
                    var relVel = vb - va;

                    double velAlongNormal = Vector2D.Dot(relVel, contact.Normal);

                    // Don't resolve if separating
                    if (velAlongNormal > 0)
                        continue;

                    // Calculate effective mass
                    double raCrossN = Vector2D.Cross(ra, contact.Normal);
                    double rbCrossN = Vector2D.Cross(rb, contact.Normal);

                    double invMassSum = a.InverseMass + b.InverseMass +
                                        raCrossN * raCrossN * a.InverseInertia +
                                        rbCrossN * rbCrossN * b.InverseInertia;

                    if (invMassSum <= PhysicsConstants.Epsilon)
                        continue;

                    // Normal impulse
                    double e = manifold.Restitution;
                    double j = -(1 + e) * velAlongNormal / invMassSum;
                    j /= manifold.ContactCount; // Distribute among contacts

                    var impulse = contact.Normal * j;
                    a.ApplyImpulseAtPoint(-impulse, contact.Point);
                    b.ApplyImpulseAtPoint(impulse, contact.Point);

                    // Friction impulse
                    va = a.Velocity + Vector2D.Cross(a.AngularVelocity, ra);
                    vb = b.Velocity + Vector2D.Cross(b.AngularVelocity, rb);
                    relVel = vb - va;

                    var tangent = relVel - contact.Normal * Vector2D.Dot(relVel, contact.Normal);
                    if (tangent.MagnitudeSquared > PhysicsConstants.Epsilon)
                    {
                        tangent = tangent.Normalized;

                        double raCrossT = Vector2D.Cross(ra, tangent);
                        double rbCrossT = Vector2D.Cross(rb, tangent);

                        double invMassSumT = a.InverseMass + b.InverseMass +
                                             raCrossT * raCrossT * a.InverseInertia +
                                             rbCrossT * rbCrossT * b.InverseInertia;

                        double jt = -Vector2D.Dot(relVel, tangent) / invMassSumT;
                        jt /= manifold.ContactCount;

                        // Coulomb friction: clamp to friction cone
                        double maxFriction = manifold.Friction * Math.Abs(j);
                        jt = Math.Clamp(jt, -maxFriction, maxFriction);

                        var frictionImpulse = tangent * jt;
                        a.ApplyImpulseAtPoint(-frictionImpulse, contact.Point);
                        b.ApplyImpulseAtPoint(frictionImpulse, contact.Point);
                    }
                }
            }
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
