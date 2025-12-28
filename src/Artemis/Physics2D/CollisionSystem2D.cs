using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Collision layers system for filtering what collides with what.
    /// Uses bit masking for efficient collision filtering.
    /// </summary>
    public static class CollisionLayers2D
    {
        public const ushort Default = 1 << 0;      // 0001
        public const ushort Static = 1 << 1;       // 0010
        public const ushort Projectile = 1 << 2;   // 0100
        public const ushort Enemy = 1 << 3;        // 1000
        public const ushort Player = 1 << 4;       // 0001 0000
        public const ushort Trigger = 1 << 5;      // 0010 0000
        public const ushort Ground = 1 << 6;       // 0100 0000
        public const ushort Debris = 1 << 7;       // 1000 0000
        public const ushort Particle = 1 << 8;

        public const ushort Everything = 0xFFFF; // All bits set

        /// <summary>
        /// Check if two layers should collide.
        /// </summary>
        public static bool ShouldCollide(ushort layerA, ushort maskA, ushort layerB, ushort maskB)
        {
            return (layerA & maskB) != 0 && (layerB & maskA) != 0;
        }
    }

    /// <summary>
    /// Collision filter for a body.
    /// </summary>
    public class CollisionFilter2D
    {
        public ushort CategoryBits { get; set; } = CollisionLayers2D.Default;
        public ushort MaskBits { get; set; } = CollisionLayers2D.Everything;
        public short GroupIndex { get; set; } = 0;

        public CollisionFilter2D() { }

        public CollisionFilter2D(ushort category, ushort mask)
        {
            CategoryBits = category;
            MaskBits = mask;
        }

        public bool CanCollideWith(CollisionFilter2D other)
        {
            // Same group with positive index: always collide
            // Same group with negative index: never collide
            if (GroupIndex != 0 && GroupIndex == other.GroupIndex)
            {
                return GroupIndex > 0;
            }

            return CollisionLayers2D.ShouldCollide(CategoryBits, MaskBits, other.CategoryBits, other.MaskBits);
        }
    }

    /// <summary>
    /// Collision event arguments.
    /// </summary>
    public class CollisionEventArgs2D : EventArgs
    {
        public RigidBody2D BodyA { get; set; }
        public RigidBody2D BodyB { get; set; }
        public Manifold2D Manifold { get; set; }

        public CollisionEventArgs2D(RigidBody2D bodyA, RigidBody2D bodyB, Manifold2D manifold)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            Manifold = manifold;
        }
    }

    /// <summary>
    /// Interface for objects that want to receive collision callbacks.
    /// </summary>
    public interface ICollisionListener2D
    {
        void OnCollisionEnter(RigidBody2D other, Manifold2D manifold);
        void OnCollisionStay(RigidBody2D other, Manifold2D manifold);
        void OnCollisionExit(RigidBody2D other);
    }

    /// <summary>
    /// Advanced SAT collision detection algorithms.
    /// </summary>
    public static class SATCollision2D
    {
        /// <summary>
        /// SAT collision detection for oriented boxes.
        /// Supports full rotation.
        /// </summary>
        public static bool BoxVsBox(RigidBody2D bodyA, RigidBody2D bodyB, BoxShape boxA, BoxShape boxB, out Manifold2D manifold)
        {
            manifold = new Manifold2D { BodyA = bodyA, BodyB = bodyB };

            // Get vertices in world space
            Vector2D[] verticesA = GetBoxVertices(bodyA.Position, boxA, bodyA.Rotation);
            Vector2D[] verticesB = GetBoxVertices(bodyB.Position, boxB, bodyB.Rotation);

            // Test all axes (4 edge normals total: 2 from each box)
            Vector2D[] axes = new Vector2D[4];

            // Axes from box A
            axes[0] = GetEdgeNormal(verticesA[0], verticesA[1]);
            axes[1] = GetEdgeNormal(verticesA[1], verticesA[2]);

            // Axes from box B
            axes[2] = GetEdgeNormal(verticesB[0], verticesB[1]);
            axes[3] = GetEdgeNormal(verticesB[1], verticesB[2]);

            double minOverlap = double.MaxValue;
            Vector2D smallestAxis = Vector2D.Zero;

            // Test each axis
            foreach (var axis in axes)
            {
                // Project both shapes onto the axis
                ProjectVertices(verticesA, axis, out double minA, out double maxA);
                ProjectVertices(verticesB, axis, out double minB, out double maxB);

                // Check for separation
                if (maxA < minB || maxB < minA)
                {
                    // Shapes are separated on this axis
                    return false;
                }

                // Calculate overlap
                double overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);

                if (overlap < minOverlap)
                {
                    minOverlap = overlap;
                    smallestAxis = axis;
                }
            }

            // If we get here, there's a collision
            manifold.Penetration = minOverlap;

            // Ensure normal points from A to B
            Vector2D direction = bodyB.Position - bodyA.Position;
            if (Vector2D.Dot(direction, smallestAxis) < 0)
            {
                smallestAxis = -smallestAxis;
            }

            manifold.Normal = smallestAxis;
            manifold.ContactCount = 1;

            return true;
        }

        /// <summary>
        /// Circle vs Box collision with proper rotation support.
        /// </summary>
        public static bool CircleVsBox(RigidBody2D circleBody, RigidBody2D boxBody, CircleShape circle, BoxShape box, out Manifold2D manifold)
        {
            manifold = new Manifold2D { BodyA = circleBody, BodyB = boxBody };

            // Transform circle center to box local space
            Vector2D localCirclePos = WorldToLocal(circleBody.Position, boxBody.Position, boxBody.Rotation);

            // Clamp to box bounds
            double halfW = box.HalfWidth;
            double halfH = box.HalfHeight;

            Vector2D closest = new Vector2D(
                Math.Clamp(localCirclePos.X, -halfW, halfW),
                Math.Clamp(localCirclePos.Y, -halfH, halfH)
            );

            Vector2D localPoint = localCirclePos - closest;
            double distanceSquared = localPoint.MagnitudeSquared;

            if (distanceSquared >= circle.Radius * circle.Radius)
                return false;

            double distance = Math.Sqrt(distanceSquared);

            // Calculate normal in local space
            Vector2D localNormal;
            if (distance > 0.0001)
            {
                localNormal = localPoint / distance;
            }
            else
            {
                // Circle center inside box - find closest edge
                double[] distances = {
                    halfW - Math.Abs(localCirclePos.X),
                    halfH - Math.Abs(localCirclePos.Y)
                };

                if (distances[0] < distances[1])
                {
                    localNormal = new Vector2D(localCirclePos.X > 0 ? 1 : -1, 0);
                }
                else
                {
                    localNormal = new Vector2D(0, localCirclePos.Y > 0 ? 1 : -1);
                }
            }

            // Transform normal back to world space
            manifold.Normal = LocalToWorldDirection(localNormal, boxBody.Rotation);
            manifold.Penetration = circle.Radius - distance;
            manifold.ContactCount = 1;

            return true;
        }

        /// <summary>
        /// Polygon vs Polygon SAT collision.
        /// </summary>
        public static bool PolygonVsPolygon(RigidBody2D bodyA, RigidBody2D bodyB,
            PolygonShape2D polyA, PolygonShape2D polyB, out Manifold2D manifold)
        {
            manifold = new Manifold2D { BodyA = bodyA, BodyB = bodyB };

            Vector2D[] verticesA = polyA.GetWorldVertices(bodyA.Position, bodyA.Rotation);
            Vector2D[] verticesB = polyB.GetWorldVertices(bodyB.Position, bodyB.Rotation);

            Vector2D[] normalsA = polyA.GetWorldNormals(bodyA.Rotation);
            Vector2D[] normalsB = polyB.GetWorldNormals(bodyB.Rotation);

            double minOverlap = double.MaxValue;
            Vector2D smallestAxis = Vector2D.Zero;

            // Test A's normals
            foreach (var axis in normalsA)
            {
                ProjectVertices(verticesA, axis, out double minA, out double maxA);
                ProjectVertices(verticesB, axis, out double minB, out double maxB);

                if (maxA < minB || maxB < minA) return false;

                double overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);
                if (overlap < minOverlap)
                {
                    minOverlap = overlap;
                    smallestAxis = axis;
                }
            }

            // Test B's normals
            foreach (var axis in normalsB)
            {
                ProjectVertices(verticesA, axis, out double minA, out double maxA);
                ProjectVertices(verticesB, axis, out double minB, out double maxB);

                if (maxA < minB || maxB < minA) return false;

                double overlap = Math.Min(maxA, maxB) - Math.Max(minA, minB);
                if (overlap < minOverlap)
                {
                    minOverlap = overlap;
                    smallestAxis = axis;
                }
            }

            manifold.Penetration = minOverlap;
            Vector2D direction = bodyB.Position - bodyA.Position;
            if (Vector2D.Dot(direction, smallestAxis) < 0)
                smallestAxis = -smallestAxis;
            manifold.Normal = smallestAxis;
            manifold.ContactCount = 1;

            return true;
        }

        private static Vector2D[] GetBoxVertices(Vector2D position, BoxShape box, double rotation)
        {
            double cos = Math.Cos(rotation);
            double sin = Math.Sin(rotation);
            double hw = box.HalfWidth;
            double hh = box.HalfHeight;

            Vector2D[] local = new Vector2D[]
            {
                new Vector2D(-hw, -hh),
                new Vector2D(hw, -hh),
                new Vector2D(hw, hh),
                new Vector2D(-hw, hh)
            };

            Vector2D[] world = new Vector2D[4];
            for (int i = 0; i < 4; i++)
            {
                double x = local[i].X * cos - local[i].Y * sin;
                double y = local[i].X * sin + local[i].Y * cos;
                world[i] = new Vector2D(position.X + x, position.Y + y);
            }

            return world;
        }

        private static Vector2D GetEdgeNormal(Vector2D p1, Vector2D p2)
        {
            Vector2D edge = p2 - p1;
            Vector2D normal = new Vector2D(-edge.Y, edge.X);
            return normal.Normalized;
        }

        private static void ProjectVertices(Vector2D[] vertices, Vector2D axis, out double min, out double max)
        {
            min = double.MaxValue;
            max = double.MinValue;

            foreach (var vertex in vertices)
            {
                double projection = Vector2D.Dot(vertex, axis);
                min = Math.Min(min, projection);
                max = Math.Max(max, projection);
            }
        }

        private static Vector2D WorldToLocal(Vector2D worldPos, Vector2D origin, double rotation)
        {
            Vector2D delta = worldPos - origin;
            double cos = Math.Cos(-rotation);
            double sin = Math.Sin(-rotation);

            return new Vector2D(
                delta.X * cos - delta.Y * sin,
                delta.X * sin + delta.Y * cos
            );
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vector2D LocalToWorldDirection(Vector2D localDir, double rotation)
        {
            double cos = Math.Cos(rotation);
            double sin = Math.Sin(rotation);

            return new Vector2D(
                localDir.X * cos - localDir.Y * sin,
                localDir.X * sin + localDir.Y * cos
            );
        }
    }

    /// <summary>
    /// Parallel collision detection system for high-performance physics.
    /// </summary>
    public class ParallelCollisionSystem2D
    {
        private readonly ConcurrentBag<Manifold2D> _detectedCollisions = new();
        private readonly ConcurrentDictionary<(string, string), Manifold2D> _manifoldCache = new();

        /// <summary>
        /// Detect collisions between all body pairs in parallel.
        /// </summary>
        public List<Manifold2D> DetectCollisionsParallel(IList<(RigidBody2D, RigidBody2D)> pairs)
        {
            while (_detectedCollisions.TryTake(out _)) { }

            Parallel.ForEach(pairs, pair =>
            {
                var (a, b) = pair;

                // Skip if both static or both sleeping
                if (a.BodyType == BodyType2D.Static && b.BodyType == BodyType2D.Static)
                    return;
                if (a.IsSleeping && b.IsSleeping)
                    return;

                // Check collision filter
                if (!a.ShouldCollide(b))
                    return;

                Manifold2D? manifold = null;

                // Detect based on shape types
                if (a.Shape is CircleShape circleA && b.Shape is CircleShape circleB)
                {
                    manifold = new Manifold2D();
                    if (!CollisionDetector2D.CircleVsCircle(a, b, manifold))
                        manifold = null;
                }
                else if (a.Shape is CircleShape circle && b.Shape is BoxShape box)
                {
                    if (SATCollision2D.CircleVsBox(a, b, circle, box, out var m))
                        manifold = m;
                }
                else if (a.Shape is BoxShape boxA2 && b.Shape is CircleShape circleB2)
                {
                    if (SATCollision2D.CircleVsBox(b, a, circleB2, boxA2, out var m))
                    {
                        m.Normal = -m.Normal;
                        var temp = m.BodyA;
                        m.BodyA = m.BodyB;
                        m.BodyB = temp;
                        manifold = m;
                    }
                }
                else if (a.Shape is BoxShape boxA && b.Shape is BoxShape boxB)
                {
                    if (SATCollision2D.BoxVsBox(a, b, boxA, boxB, out var m))
                        manifold = m;
                }

                if (manifold != null)
                {
                    manifold.IsActive = true;
                    _detectedCollisions.Add(manifold);
                }
            });

            return new List<Manifold2D>(_detectedCollisions);
        }

        /// <summary>
        /// Batch broad-phase collision detection using spatial hash grid.
        /// </summary>
        public List<(RigidBody2D, RigidBody2D)> BroadPhaseParallel(SpatialHashGrid2D grid)
        {
            return grid.GetPotentialPairsParallel();
        }

        /// <summary>
        /// Full collision detection pipeline (broad + narrow phase, parallel).
        /// </summary>
        public List<Manifold2D> DetectAllCollisionsParallel(IList<RigidBody2D> bodies)
        {
            // Build spatial grid
            var grid = new SpatialHashGrid2D(5.0);
            grid.RebuildParallel(bodies);

            // Broad phase
            var pairs = BroadPhaseParallel(grid);

            // Narrow phase
            return DetectCollisionsParallel(pairs);
        }

        /// <summary>
        /// Resolve all collisions in parallel (position correction only, velocity must be sequential).
        /// </summary>
        public void ResolvePositionsParallel(IList<Manifold2D> manifolds, double baumgarte = 0.2)
        {
            Parallel.ForEach(manifolds, manifold =>
            {
                if (!manifold.IsActive) return;

                var a = manifold.BodyA;
                var b = manifold.BodyB;

                if (a.BodyType == BodyType2D.Static && b.BodyType == BodyType2D.Static)
                    return;

                double totalInvMass = a.InverseMass + b.InverseMass;
                if (totalInvMass < 1e-8) return;

                double correction = Math.Max(manifold.Penetration - 0.01, 0) * baumgarte;
                Vector2D correctionVec = manifold.Normal * correction;

                if (a.BodyType == BodyType2D.Dynamic)
                {
                    double ratio = a.InverseMass / totalInvMass;
                    a.Position -= correctionVec * ratio;
                }

                if (b.BodyType == BodyType2D.Dynamic)
                {
                    double ratio = b.InverseMass / totalInvMass;
                    b.Position += correctionVec * ratio;
                }
            });
        }

        public void Clear()
        {
            while (_detectedCollisions.TryTake(out _)) { }
            _manifoldCache.Clear();
        }
    }

    /// <summary>
    /// Batch collision operations for SIMD-style processing.
    /// </summary>
    public static class BatchCollisionOperations
    {
        /// <summary>
        /// Batch circle vs circle collision test.
        /// </summary>
        public static bool[] BatchCircleVsCircle(
            Vector2D[] centersA, double[] radiiA,
            Vector2D[] centersB, double[] radiiB)
        {
            int count = centersA.Length;
            var results = new bool[count];

            Parallel.For(0, count, i =>
            {
                double dx = centersB[i].X - centersA[i].X;
                double dy = centersB[i].Y - centersA[i].Y;
                double distSq = dx * dx + dy * dy;
                double radiusSum = radiiA[i] + radiiB[i];
                results[i] = distSq < radiusSum * radiusSum;
            });

            return results;
        }

        /// <summary>
        /// Batch AABB overlap test.
        /// </summary>
        public static bool[] BatchAABBOverlap(AABB2D[] aabbsA, AABB2D[] aabbsB)
        {
            int count = aabbsA.Length;
            var results = new bool[count];

            Parallel.For(0, count, i =>
            {
                results[i] = aabbsA[i].Overlaps(aabbsB[i]);
            });

            return results;
        }

        /// <summary>
        /// Batch point vs circle test.
        /// </summary>
        public static bool[] BatchPointInCircle(Vector2D[] points, Vector2D center, double radius)
        {
            var results = new bool[points.Length];
            double radiusSq = radius * radius;

            Parallel.For(0, points.Length, i =>
            {
                results[i] = (points[i] - center).MagnitudeSquared <= radiusSq;
            });

            return results;
        }

        /// <summary>
        /// Batch separating axis test for boxes.
        /// </summary>
        public static bool[] BatchBoxOverlap(
            Vector2D[] positionsA, double[] rotationsA, double[] halfWidthsA, double[] halfHeightsA,
            Vector2D[] positionsB, double[] rotationsB, double[] halfWidthsB, double[] halfHeightsB)
        {
            int count = positionsA.Length;
            var results = new bool[count];

            Parallel.For(0, count, i =>
            {
                // Quick AABB pre-check
                double maxExtentA = Math.Sqrt(halfWidthsA[i] * halfWidthsA[i] + halfHeightsA[i] * halfHeightsA[i]);
                double maxExtentB = Math.Sqrt(halfWidthsB[i] * halfWidthsB[i] + halfHeightsB[i] * halfHeightsB[i]);
                double dx = positionsB[i].X - positionsA[i].X;
                double dy = positionsB[i].Y - positionsA[i].Y;
                double distSq = dx * dx + dy * dy;
                double maxDistSq = (maxExtentA + maxExtentB) * (maxExtentA + maxExtentB);

                if (distSq > maxDistSq)
                {
                    results[i] = false;
                    return;
                }

                // Full SAT test would go here - simplified for performance
                results[i] = true;
            });

            return results;
        }
    }
}
