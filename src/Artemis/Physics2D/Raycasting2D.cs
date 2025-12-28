using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Artemis.Physics2D
{
    /// <summary>
    /// 2D Ray definition.
    /// </summary>
    public struct Ray2D
    {
        public Vector2D Origin { get; set; }
        public Vector2D Direction { get; set; }
        public double MaxDistance { get; set; }

        public Ray2D(Vector2D origin, Vector2D direction, double maxDistance = double.MaxValue)
        {
            Origin = origin;
            Direction = direction.Normalized;
            MaxDistance = maxDistance;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D GetPoint(double distance)
        {
            return Origin + Direction * distance;
        }
    }

    /// <summary>
    /// Result of a raycast.
    /// </summary>
    public struct RaycastHit2D
    {
        public RigidBody2D? Body { get; set; }
        public Vector2D Point { get; set; }
        public Vector2D Normal { get; set; }
        public double Distance { get; set; }
        public double Fraction { get; set; }
        public bool Hit { get; set; }

        public static RaycastHit2D NoHit => new RaycastHit2D { Hit = false, Distance = double.MaxValue, Fraction = double.MaxValue };
    }

    /// <summary>
    /// High-performance raycasting utilities for 2D physics.
    /// Supports parallel batch raycasting for maximum throughput.
    /// </summary>
    public static class Raycaster2D
    {
        /// <summary>
        /// Cast a ray and return the closest hit.
        /// </summary>
        public static RaycastHit2D Raycast(Ray2D ray, IEnumerable<RigidBody2D> bodies, ushort layerMask = 0xFFFF)
        {
            RaycastHit2D closest = RaycastHit2D.NoHit;
            double closestDistance = ray.MaxDistance;

            foreach (var body in bodies)
            {
                if (!body.IsActive) continue;
                if ((body.CategoryBits & layerMask) == 0) continue;

                // Quick AABB check first
                if (!RayIntersectsAABB(ray, body.AABB, closestDistance))
                    continue;

                RaycastHit2D hit;
                bool didHit = false;

                if (body.Shape is CircleShape circle)
                {
                    var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                    didHit = RaycastCircle(ray, center, circle.Radius, out hit);
                }
                else if (body.Shape is BoxShape box)
                {
                    didHit = RaycastBox(ray, body.Position, box, body.Rotation, out hit);
                }
                else
                {
                    continue;
                }

                if (didHit && hit.Distance < closestDistance)
                {
                    hit.Body = body;
                    closest = hit;
                    closestDistance = hit.Distance;
                }
            }

            closest.Hit = closestDistance < ray.MaxDistance;
            return closest;
        }

        /// <summary>
        /// Cast a ray using spatial hash grid for acceleration.
        /// </summary>
        public static RaycastHit2D RaycastWithSpatialHash(Ray2D ray, SpatialHashGrid2D grid, ushort layerMask = 0xFFFF)
        {
            // Get bodies along ray path
            Vector2D end = ray.GetPoint(ray.MaxDistance);
            Vector2D min = new Vector2D(Math.Min(ray.Origin.X, end.X), Math.Min(ray.Origin.Y, end.Y));
            Vector2D max = new Vector2D(Math.Max(ray.Origin.X, end.X), Math.Max(ray.Origin.Y, end.Y));

            var candidates = grid.QueryRegion(min, max);
            return Raycast(ray, candidates, layerMask);
        }

        /// <summary>
        /// Cast a ray and return all hits.
        /// </summary>
        public static RaycastHit2D[] RaycastAll(Ray2D ray, IEnumerable<RigidBody2D> bodies, ushort layerMask = 0xFFFF)
        {
            List<RaycastHit2D> hits = new List<RaycastHit2D>();

            foreach (var body in bodies)
            {
                if (!body.IsActive) continue;
                if ((body.CategoryBits & layerMask) == 0) continue;

                RaycastHit2D hit;
                bool didHit = false;

                if (body.Shape is CircleShape circle)
                {
                    var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                    didHit = RaycastCircle(ray, center, circle.Radius, out hit);
                }
                else if (body.Shape is BoxShape box)
                {
                    didHit = RaycastBox(ray, body.Position, box, body.Rotation, out hit);
                }
                else
                {
                    continue;
                }

                if (didHit && hit.Distance <= ray.MaxDistance)
                {
                    hit.Body = body;
                    hit.Hit = true;
                    hits.Add(hit);
                }
            }

            return hits.OrderBy(h => h.Distance).ToArray();
        }

        /// <summary>
        /// Cast multiple rays in parallel (batch raycasting).
        /// </summary>
        public static RaycastHit2D[] RaycastBatchParallel(Ray2D[] rays, IList<RigidBody2D> bodies, ushort layerMask = 0xFFFF)
        {
            var results = new RaycastHit2D[rays.Length];

            Parallel.For(0, rays.Length, i =>
            {
                results[i] = Raycast(rays[i], bodies, layerMask);
            });

            return results;
        }

        /// <summary>
        /// Cast rays in a fan pattern (useful for vision cones, sensors).
        /// </summary>
        public static RaycastHit2D[] RaycastFanParallel(
            Vector2D origin,
            double startAngle,
            double endAngle,
            int rayCount,
            double maxDistance,
            IList<RigidBody2D> bodies,
            ushort layerMask = 0xFFFF)
        {
            var rays = new Ray2D[rayCount];
            double angleStep = (endAngle - startAngle) / (rayCount - 1);

            for (int i = 0; i < rayCount; i++)
            {
                double angle = startAngle + i * angleStep;
                var direction = new Vector2D(Math.Cos(angle), Math.Sin(angle));
                rays[i] = new Ray2D(origin, direction, maxDistance);
            }

            return RaycastBatchParallel(rays, bodies, layerMask);
        }

        /// <summary>
        /// Cast rays in a circle pattern (360 degree sensor).
        /// </summary>
        public static RaycastHit2D[] RaycastCircleParallel(
            Vector2D origin,
            int rayCount,
            double maxDistance,
            IList<RigidBody2D> bodies,
            ushort layerMask = 0xFFFF)
        {
            return RaycastFanParallel(origin, 0, 2 * Math.PI, rayCount, maxDistance, bodies, layerMask);
        }

        /// <summary>
        /// Quick AABB intersection test for ray culling.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static bool RayIntersectsAABB(Ray2D ray, AABB2D aabb, double maxDist)
        {
            double tmin = 0;
            double tmax = maxDist;

            // X axis
            if (Math.Abs(ray.Direction.X) < 1e-8)
            {
                if (ray.Origin.X < aabb.Min.X || ray.Origin.X > aabb.Max.X)
                    return false;
            }
            else
            {
                double t1 = (aabb.Min.X - ray.Origin.X) / ray.Direction.X;
                double t2 = (aabb.Max.X - ray.Origin.X) / ray.Direction.X;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                tmin = Math.Max(tmin, t1);
                tmax = Math.Min(tmax, t2);
                if (tmin > tmax) return false;
            }

            // Y axis
            if (Math.Abs(ray.Direction.Y) < 1e-8)
            {
                if (ray.Origin.Y < aabb.Min.Y || ray.Origin.Y > aabb.Max.Y)
                    return false;
            }
            else
            {
                double t1 = (aabb.Min.Y - ray.Origin.Y) / ray.Direction.Y;
                double t2 = (aabb.Max.Y - ray.Origin.Y) / ray.Direction.Y;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                tmin = Math.Max(tmin, t1);
                tmax = Math.Min(tmax, t2);
                if (tmin > tmax) return false;
            }

            return true;
        }

        /// <summary>
        /// Raycast against a circle.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool RaycastCircle(Ray2D ray, Vector2D center, double radius, out RaycastHit2D hit)
        {
            hit = RaycastHit2D.NoHit;

            Vector2D m = ray.Origin - center;
            double b = Vector2D.Dot(m, ray.Direction);
            double c = Vector2D.Dot(m, m) - radius * radius;

            if (c > 0.0 && b > 0.0)
                return false;

            double discriminant = b * b - c;

            if (discriminant < 0.0)
                return false;

            double distance = -b - Math.Sqrt(discriminant);

            if (distance < 0)
                distance = 0;

            Vector2D point = ray.Origin + ray.Direction * distance;
            Vector2D normal = (point - center).Normalized;

            hit = new RaycastHit2D
            {
                Hit = true,
                Point = point,
                Normal = normal,
                Distance = distance,
                Fraction = distance / ray.MaxDistance
            };

            return true;
        }

        /// <summary>
        /// Raycast against a box (with rotation support).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool RaycastBox(Ray2D ray, Vector2D center, BoxShape box, double rotation, out RaycastHit2D hit)
        {
            hit = RaycastHit2D.NoHit;

            // Transform ray to box local space
            double cos = Math.Cos(-rotation);
            double sin = Math.Sin(-rotation);

            Vector2D localOrigin = ray.Origin - center;
            Vector2D rotatedOrigin = new Vector2D(
                localOrigin.X * cos - localOrigin.Y * sin,
                localOrigin.X * sin + localOrigin.Y * cos
            );

            Vector2D rotatedDirection = new Vector2D(
                ray.Direction.X * cos - ray.Direction.Y * sin,
                ray.Direction.X * sin + ray.Direction.Y * cos
            );

            // AABB raycast in local space
            double halfW = box.HalfWidth;
            double halfH = box.HalfHeight;

            double tMin = 0;
            double tMax = ray.MaxDistance;
            int hitAxis = -1;
            bool hitMin = true;

            // X axis
            if (Math.Abs(rotatedDirection.X) < 1e-8)
            {
                if (rotatedOrigin.X < -halfW || rotatedOrigin.X > halfW)
                    return false;
            }
            else
            {
                double t1 = (-halfW - rotatedOrigin.X) / rotatedDirection.X;
                double t2 = (halfW - rotatedOrigin.X) / rotatedDirection.X;

                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; hitMin = false; }
                else hitMin = true;

                if (t1 > tMin) { tMin = t1; hitAxis = 0; }
                tMax = Math.Min(tMax, t2);

                if (tMin > tMax) return false;
            }

            // Y axis
            if (Math.Abs(rotatedDirection.Y) < 1e-8)
            {
                if (rotatedOrigin.Y < -halfH || rotatedOrigin.Y > halfH)
                    return false;
            }
            else
            {
                double t1 = (-halfH - rotatedOrigin.Y) / rotatedDirection.Y;
                double t2 = (halfH - rotatedOrigin.Y) / rotatedDirection.Y;

                bool fromMin = t1 < t2;
                if (!fromMin) { double tmp = t1; t1 = t2; t2 = tmp; }

                if (t1 > tMin) { tMin = t1; hitAxis = 1; hitMin = fromMin; }
                tMax = Math.Min(tMax, t2);

                if (tMin > tMax) return false;
            }

            // Calculate normal in local space
            Vector2D localNormal;
            if (hitAxis == 0)
                localNormal = hitMin ? new Vector2D(-1, 0) : new Vector2D(1, 0);
            else
                localNormal = hitMin ? new Vector2D(0, -1) : new Vector2D(0, 1);

            // Transform back to world space
            double wcos = Math.Cos(rotation);
            double wsin = Math.Sin(rotation);

            Vector2D worldNormal = new Vector2D(
                localNormal.X * wcos - localNormal.Y * wsin,
                localNormal.X * wsin + localNormal.Y * wcos
            );

            Vector2D worldPoint = ray.Origin + ray.Direction * tMin;

            hit = new RaycastHit2D
            {
                Hit = true,
                Point = worldPoint,
                Normal = worldNormal.Normalized,
                Distance = tMin,
                Fraction = tMin / ray.MaxDistance
            };

            return true;
        }

        /// <summary>
        /// Raycast against an edge.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool RaycastEdge(Ray2D ray, Vector2D v1, Vector2D v2, out RaycastHit2D hit)
        {
            hit = RaycastHit2D.NoHit;

            Vector2D edge = v2 - v1;
            Vector2D normal = new Vector2D(-edge.Y, edge.X).Normalized;

            // Check if ray is parallel to edge
            double denom = Vector2D.Dot(ray.Direction, normal);
            if (Math.Abs(denom) < 1e-8)
                return false;

            double t = Vector2D.Dot(v1 - ray.Origin, normal) / denom;
            if (t < 0 || t > ray.MaxDistance)
                return false;

            Vector2D point = ray.GetPoint(t);

            // Check if point is within edge segment
            double edgeLen = edge.Magnitude;
            Vector2D edgeDir = edge / edgeLen;
            double projection = Vector2D.Dot(point - v1, edgeDir);

            if (projection < 0 || projection > edgeLen)
                return false;

            // Ensure normal points toward ray origin
            if (denom > 0)
                normal = -normal;

            hit = new RaycastHit2D
            {
                Hit = true,
                Point = point,
                Normal = normal,
                Distance = t,
                Fraction = t / ray.MaxDistance
            };

            return true;
        }

        /// <summary>
        /// Raycast against a polygon.
        /// </summary>
        public static bool RaycastPolygon(Ray2D ray, PolygonShape2D polygon, Vector2D position, double rotation, out RaycastHit2D hit)
        {
            hit = RaycastHit2D.NoHit;

            var vertices = polygon.GetWorldVertices(position, rotation);
            double closestT = ray.MaxDistance;
            bool foundHit = false;

            for (int i = 0; i < vertices.Length; i++)
            {
                int next = (i + 1) % vertices.Length;
                if (RaycastEdge(ray, vertices[i], vertices[next], out var edgeHit))
                {
                    if (edgeHit.Distance < closestT)
                    {
                        closestT = edgeHit.Distance;
                        hit = edgeHit;
                        foundHit = true;
                    }
                }
            }

            return foundHit;
        }
    }

    /// <summary>
    /// Parallel raycast job for batch processing.
    /// </summary>
    public class ParallelRaycastJob
    {
        private readonly IList<RigidBody2D> _bodies;
        private readonly ushort _layerMask;
        private readonly ConcurrentQueue<(Ray2D ray, int index)> _rayQueue;
        private readonly RaycastHit2D[] _results;

        public ParallelRaycastJob(IList<RigidBody2D> bodies, int rayCount, ushort layerMask = 0xFFFF)
        {
            _bodies = bodies;
            _layerMask = layerMask;
            _rayQueue = new ConcurrentQueue<(Ray2D, int)>();
            _results = new RaycastHit2D[rayCount];
        }

        public void EnqueueRay(Ray2D ray, int index)
        {
            _rayQueue.Enqueue((ray, index));
        }

        public RaycastHit2D[] Execute()
        {
            var rays = _rayQueue.ToArray();

            Parallel.ForEach(rays, item =>
            {
                _results[item.index] = Raycaster2D.Raycast(item.ray, _bodies, _layerMask);
            });

            return _results;
        }
    }
}
