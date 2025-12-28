using System;
using System.Collections.Generic;
using System.Linq;

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
    /// Raycasting utilities for 2D physics.
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
        /// Raycast against a circle.
        /// </summary>
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

            Vector2D min = new Vector2D(-halfW, -halfH);
            Vector2D max = new Vector2D(halfW, halfH);

            double tMin = 0;
            double tMax = ray.MaxDistance;
            int hitAxis = -1;

            for (int i = 0; i < 2; i++)
            {
                double origin = i == 0 ? rotatedOrigin.X : rotatedOrigin.Y;
                double dir = i == 0 ? rotatedDirection.X : rotatedDirection.Y;
                double minVal = i == 0 ? min.X : min.Y;
                double maxVal = i == 0 ? max.X : max.Y;

                if (Math.Abs(dir) < 0.0001)
                {
                    if (origin < minVal || origin > maxVal)
                        return false;
                }
                else
                {
                    double t1 = (minVal - origin) / dir;
                    double t2 = (maxVal - origin) / dir;

                    bool fromMin = t1 < t2;
                    if (!fromMin)
                    {
                        double temp = t1;
                        t1 = t2;
                        t2 = temp;
                    }

                    if (t1 > tMin)
                    {
                        tMin = t1;
                        hitAxis = i;
                    }
                    tMax = Math.Min(tMax, t2);

                    if (tMin > tMax)
                        return false;
                }
            }

            Vector2D localPoint = rotatedOrigin + rotatedDirection * tMin;

            // Calculate normal in local space
            Vector2D localNormal;
            if (hitAxis == 0)
            {
                localNormal = rotatedDirection.X > 0 ? new Vector2D(-1, 0) : new Vector2D(1, 0);
            }
            else
            {
                localNormal = rotatedDirection.Y > 0 ? new Vector2D(0, -1) : new Vector2D(0, 1);
            }

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
        public static bool RaycastEdge(Ray2D ray, Vector2D v1, Vector2D v2, out RaycastHit2D hit)
        {
            hit = RaycastHit2D.NoHit;

            Vector2D edge = v2 - v1;
            Vector2D normal = new Vector2D(-edge.Y, edge.X).Normalized;

            // Check if ray is parallel to edge
            double denom = Vector2D.Dot(ray.Direction, normal);
            if (Math.Abs(denom) < 0.0001)
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
    }
}
