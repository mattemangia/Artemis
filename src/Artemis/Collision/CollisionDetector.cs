using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Collision
{
    /// <summary>
    /// Provides collision detection algorithms for various shape pairs.
    /// </summary>
    public static class CollisionDetector
    {
        #region Main Detection Method

        /// <summary>
        /// Detects collision between two physics bodies.
        /// </summary>
        public static CollisionInfo Detect(IPhysicsBody bodyA, IPhysicsBody bodyB)
        {
            // Broad phase: AABB test
            if (!bodyA.BoundingBox.Intersects(bodyB.BoundingBox))
                return CollisionInfo.NoCollision;

            // Get rigid bodies for shape info
            if (bodyA is RigidBody rbA && bodyB is RigidBody rbB)
            {
                return DetectRigidBodies(rbA, rbB);
            }

            return CollisionInfo.NoCollision;
        }

        private static CollisionInfo DetectRigidBodies(RigidBody bodyA, RigidBody bodyB)
        {
            return (bodyA.ShapeType, bodyB.ShapeType) switch
            {
                (CollisionShapeType.Sphere, CollisionShapeType.Sphere) => SphereSphere(bodyA, bodyB),
                (CollisionShapeType.Sphere, CollisionShapeType.Box) => SphereBox(bodyA, bodyB),
                (CollisionShapeType.Box, CollisionShapeType.Sphere) => SphereBox(bodyB, bodyA).Swapped(),
                (CollisionShapeType.Box, CollisionShapeType.Box) => BoxBox(bodyA, bodyB),
                _ => CollisionInfo.NoCollision
            };
        }

        #endregion

        #region Sphere-Sphere

        /// <summary>
        /// Detects collision between two spheres.
        /// </summary>
        public static CollisionInfo SphereSphere(RigidBody sphereA, RigidBody sphereB)
        {
            var delta = sphereB.Position - sphereA.Position;
            double distanceSquared = delta.MagnitudeSquared;
            double radiusSum = sphereA.Radius + sphereB.Radius;

            if (distanceSquared >= radiusSum * radiusSum)
                return CollisionInfo.NoCollision;

            double distance = Math.Sqrt(distanceSquared);

            Vector3D normal;
            if (distance < PhysicsConstants.Epsilon)
            {
                // Spheres are at the same position, use arbitrary normal
                normal = Vector3D.Up;
                distance = 0;
            }
            else
            {
                normal = delta / distance;
            }

            double penetration = radiusSum - distance;
            var contactPoint = sphereA.Position + normal * (sphereA.Radius - penetration * 0.5);

            return CollisionInfo.Create(sphereA, sphereB, normal, penetration, contactPoint);
        }

        #endregion

        #region Sphere-Box

        /// <summary>
        /// Detects collision between a sphere and a box.
        /// </summary>
        public static CollisionInfo SphereBox(RigidBody sphere, RigidBody box)
        {
            // Transform sphere center to box's local space
            var localCenter = box.Transform.InverseTransformPoint(sphere.Position);

            // Find closest point on box to sphere center
            var closestLocal = new Vector3D(
                Math.Clamp(localCenter.X, -box.HalfExtents.X, box.HalfExtents.X),
                Math.Clamp(localCenter.Y, -box.HalfExtents.Y, box.HalfExtents.Y),
                Math.Clamp(localCenter.Z, -box.HalfExtents.Z, box.HalfExtents.Z)
            );

            // Check if sphere center is inside the box
            bool inside = (localCenter - closestLocal).MagnitudeSquared < PhysicsConstants.Epsilon;

            Vector3D normal;
            double penetration;
            Vector3D contactPoint;

            if (inside)
            {
                // Sphere center is inside box, find nearest face
                var distances = new double[]
                {
                    box.HalfExtents.X - Math.Abs(localCenter.X),
                    box.HalfExtents.Y - Math.Abs(localCenter.Y),
                    box.HalfExtents.Z - Math.Abs(localCenter.Z)
                };

                int minAxis = 0;
                double minDist = distances[0];
                for (int i = 1; i < 3; i++)
                {
                    if (distances[i] < minDist)
                    {
                        minDist = distances[i];
                        minAxis = i;
                    }
                }

                // Push out along the nearest face
                var localNormal = Vector3D.Zero;
                switch (minAxis)
                {
                    case 0:
                        localNormal = new Vector3D(SignNonZero(localCenter.X), 0, 0);
                        closestLocal = new Vector3D(
                            SignNonZero(localCenter.X) * box.HalfExtents.X,
                            localCenter.Y,
                            localCenter.Z);
                        break;
                    case 1:
                        localNormal = new Vector3D(0, SignNonZero(localCenter.Y), 0);
                        closestLocal = new Vector3D(
                            localCenter.X,
                            SignNonZero(localCenter.Y) * box.HalfExtents.Y,
                            localCenter.Z);
                        break;
                    case 2:
                        localNormal = new Vector3D(0, 0, SignNonZero(localCenter.Z));
                        closestLocal = new Vector3D(
                            localCenter.X,
                            localCenter.Y,
                            SignNonZero(localCenter.Z) * box.HalfExtents.Z);
                        break;
                }

                // The collision normal points from A (the sphere) to B (the box),
                // i.e. opposite to the direction that pushes the sphere out.
                normal = -box.Transform.TransformDirection(localNormal);
                penetration = sphere.Radius + minDist;
                contactPoint = box.Transform.TransformPoint(closestLocal);
            }
            else
            {
                // Sphere center is outside box
                var delta = localCenter - closestLocal;
                double distance = delta.Magnitude;

                if (distance >= sphere.Radius)
                    return CollisionInfo.NoCollision;

                var localNormal = delta / distance;
                normal = -box.Transform.TransformDirection(localNormal);
                penetration = sphere.Radius - distance;
                contactPoint = box.Transform.TransformPoint(closestLocal);
            }

            return CollisionInfo.Create(sphere, box, normal, penetration, contactPoint);
        }

        #endregion

        #region Box-Box (SAT)

        // Face axes are preferred over edge axes unless the edge axis is clearly
        // shallower. Face contacts give a full manifold, which keeps resting and
        // stacked boxes stable.
        private const double EdgeAxisRelativeTolerance = 0.95;
        private const double EdgeAxisAbsoluteTolerance = 0.01;

        /// <summary>
        /// Detects collision between two boxes using the Separating Axis Theorem and
        /// builds a contact manifold (up to four points for face contacts, one point
        /// for edge-edge contacts).
        /// </summary>
        public static CollisionInfo BoxBox(RigidBody boxA, RigidBody boxB)
        {
            var rotA = Matrix3x3.FromQuaternion(boxA.Rotation);
            var rotB = Matrix3x3.FromQuaternion(boxB.Rotation);

            var axesA = new Vector3D[]
            {
                new(rotA.M00, rotA.M10, rotA.M20),
                new(rotA.M01, rotA.M11, rotA.M21),
                new(rotA.M02, rotA.M12, rotA.M22)
            };

            var axesB = new Vector3D[]
            {
                new(rotB.M00, rotB.M10, rotB.M20),
                new(rotB.M01, rotB.M11, rotB.M21),
                new(rotB.M02, rotB.M12, rotB.M22)
            };

            var halfA = boxA.HalfExtents;
            var halfB = boxB.HalfExtents;
            var t = boxB.Position - boxA.Position;

            double facePenetration = double.MaxValue;
            Vector3D faceAxis = Vector3D.Zero;
            int faceAxisType = -1;

            for (int i = 0; i < 3; i++)
            {
                if (!TestAxis(axesA[i], halfA, halfB, axesA, axesB, t,
                    ref facePenetration, ref faceAxis, ref faceAxisType, i))
                    return CollisionInfo.NoCollision;
            }

            for (int i = 0; i < 3; i++)
            {
                if (!TestAxis(axesB[i], halfA, halfB, axesA, axesB, t,
                    ref facePenetration, ref faceAxis, ref faceAxisType, i + 3))
                    return CollisionInfo.NoCollision;
            }

            double edgePenetration = double.MaxValue;
            Vector3D edgeAxis = Vector3D.Zero;
            int edgeAxisType = -1;

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    var axis = Vector3D.Cross(axesA[i], axesB[j]);
                    if (axis.MagnitudeSquared < 1e-8)
                        continue; // Parallel edges: already covered by the face axes.

                    axis = axis.Normalized;
                    if (!TestAxis(axis, halfA, halfB, axesA, axesB, t,
                        ref edgePenetration, ref edgeAxis, ref edgeAxisType, 6 + i * 3 + j))
                        return CollisionInfo.NoCollision;
                }
            }

            bool useEdge = edgeAxisType >= 0 &&
                           edgePenetration < facePenetration * EdgeAxisRelativeTolerance - EdgeAxisAbsoluteTolerance;

            if (useEdge)
            {
                var normal = Vector3D.Dot(edgeAxis, t) < 0 ? -edgeAxis : edgeAxis;
                int edgeA = (edgeAxisType - 6) / 3;
                int edgeB = (edgeAxisType - 6) % 3;

                // Supporting edges: the edge of A furthest along the normal and the
                // edge of B furthest against it.
                var pointA = boxA.Position;
                var pointB = boxB.Position;
                for (int k = 0; k < 3; k++)
                {
                    if (k != edgeA)
                        pointA += axesA[k] * (SignNonZero(Vector3D.Dot(axesA[k], normal)) * Component(halfA, k));
                    if (k != edgeB)
                        pointB -= axesB[k] * (SignNonZero(Vector3D.Dot(axesB[k], normal)) * Component(halfB, k));
                }

                ClosestPointsOnLines(pointA, axesA[edgeA], Component(halfA, edgeA),
                    pointB, axesB[edgeB], Component(halfB, edgeB), out var onA, out var onB);

                return CollisionInfo.Create(boxA, boxB, normal, edgePenetration, (onA + onB) * 0.5);
            }

            // Face contact: clip the incident face against the reference face.
            bool referenceIsA = faceAxisType < 3;
            var refBody = referenceIsA ? boxA : boxB;
            var incBody = referenceIsA ? boxB : boxA;
            var refAxes = referenceIsA ? axesA : axesB;
            var incAxes = referenceIsA ? axesB : axesA;
            var refHalf = referenceIsA ? halfA : halfB;
            var incHalf = referenceIsA ? halfB : halfA;
            int refIndex = referenceIsA ? faceAxisType : faceAxisType - 3;

            // Reference face normal, pointing towards the incident box.
            var toIncident = incBody.Position - refBody.Position;
            var refNormal = refAxes[refIndex];
            if (Vector3D.Dot(refNormal, toIncident) < 0)
                refNormal = -refNormal;

            var collisionNormal = referenceIsA ? refNormal : -refNormal;

            // Incident face: the face of the other box most opposed to the reference normal.
            int incIndex = 0;
            double mostOpposed = double.MaxValue;
            double incSign = 1;
            for (int k = 0; k < 3; k++)
            {
                double d = Vector3D.Dot(incAxes[k], refNormal);
                if (-Math.Abs(d) < mostOpposed)
                {
                    mostOpposed = -Math.Abs(d);
                    incIndex = k;
                    incSign = d > 0 ? -1 : 1;
                }
            }

            int incU = (incIndex + 1) % 3;
            int incV = (incIndex + 2) % 3;
            var incCenter = incBody.Position + incAxes[incIndex] * (incSign * Component(incHalf, incIndex));
            var eu = incAxes[incU] * Component(incHalf, incU);
            var ev = incAxes[incV] * Component(incHalf, incV);

            var polygon = new List<Vector3D>(8)
            {
                incCenter + eu + ev,
                incCenter - eu + ev,
                incCenter - eu - ev,
                incCenter + eu - ev
            };
            var scratch = new List<Vector3D>(8);

            int refU = (refIndex + 1) % 3;
            int refV = (refIndex + 2) % 3;
            ClipPolygon(polygon, scratch, refAxes[refU], Vector3D.Dot(refAxes[refU], refBody.Position) + Component(refHalf, refU));
            ClipPolygon(polygon, scratch, -refAxes[refU], -Vector3D.Dot(refAxes[refU], refBody.Position) + Component(refHalf, refU));
            ClipPolygon(polygon, scratch, refAxes[refV], Vector3D.Dot(refAxes[refV], refBody.Position) + Component(refHalf, refV));
            ClipPolygon(polygon, scratch, -refAxes[refV], -Vector3D.Dot(refAxes[refV], refBody.Position) + Component(refHalf, refV));

            double refPlane = Vector3D.Dot(refNormal, refBody.Position) + Component(refHalf, refIndex);

            var points = new List<ManifoldPoint>(polygon.Count);
            foreach (var vertex in polygon)
            {
                double depth = refPlane - Vector3D.Dot(refNormal, vertex);
                if (depth >= 0)
                    points.Add(new ManifoldPoint { Position = vertex, Penetration = depth });
            }

            if (points.Count == 0)
            {
                // Numerically degenerate overlap: fall back to a single point between
                // the two faces so the pair is still resolved.
                var fallback = incCenter + refNormal * (0.5 * facePenetration);
                return CollisionInfo.Create(boxA, boxB, collisionNormal, facePenetration, fallback);
            }

            return CollisionInfo.Create(boxA, boxB, collisionNormal, ReduceManifold(points));
        }

        /// <summary>
        /// Sutherland-Hodgman clip of a convex polygon against the half space
        /// <c>dot(normal, p) &lt;= offset</c>.
        /// </summary>
        private static void ClipPolygon(List<Vector3D> polygon, List<Vector3D> scratch, Vector3D normal, double offset)
        {
            if (polygon.Count == 0)
                return;

            scratch.Clear();
            var previous = polygon[polygon.Count - 1];
            double previousDistance = Vector3D.Dot(normal, previous) - offset;

            foreach (var current in polygon)
            {
                double currentDistance = Vector3D.Dot(normal, current) - offset;

                if (currentDistance <= 0)
                {
                    if (previousDistance > 0)
                        scratch.Add(previous + (current - previous) * (previousDistance / (previousDistance - currentDistance)));
                    scratch.Add(current);
                }
                else if (previousDistance <= 0)
                {
                    scratch.Add(previous + (current - previous) * (previousDistance / (previousDistance - currentDistance)));
                }

                previous = current;
                previousDistance = currentDistance;
            }

            polygon.Clear();
            polygon.AddRange(scratch);
        }

        /// <summary>
        /// Keeps at most four manifold points: the deepest one plus the three that
        /// span the largest area around it.
        /// </summary>
        private static ManifoldPoint[] ReduceManifold(List<ManifoldPoint> points)
        {
            const int max = ContactManifold.MaxContacts;
            if (points.Count <= max)
                return points.ToArray();

            var result = new ManifoldPoint[max];
            var used = new bool[points.Count];

            int first = 0;
            for (int i = 1; i < points.Count; i++)
                if (points[i].Penetration > points[first].Penetration) first = i;
            result[0] = points[first];
            used[first] = true;

            for (int n = 1; n < max; n++)
            {
                int best = -1;
                double bestScore = -1;
                for (int i = 0; i < points.Count; i++)
                {
                    if (used[i]) continue;

                    // Distance to the closest point already chosen.
                    double score = double.MaxValue;
                    for (int k = 0; k < n; k++)
                        score = Math.Min(score, (points[i].Position - result[k].Position).MagnitudeSquared);

                    if (score > bestScore)
                    {
                        bestScore = score;
                        best = i;
                    }
                }

                result[n] = points[best];
                used[best] = true;
            }

            return result;
        }

        /// <summary>
        /// Closest points between two segments given by centre, unit direction and half length.
        /// </summary>
        private static void ClosestPointsOnLines(
            Vector3D centerA, Vector3D dirA, double halfLengthA,
            Vector3D centerB, Vector3D dirB, double halfLengthB,
            out Vector3D onA, out Vector3D onB)
        {
            var r = centerA - centerB;
            double b = Vector3D.Dot(dirA, dirB);
            double d = Vector3D.Dot(dirA, r);
            double e = Vector3D.Dot(dirB, r);
            double denom = 1.0 - b * b;

            double sA = 0, sB;
            if (denom > 1e-9)
                sA = Math.Clamp((b * e - d) / denom, -halfLengthA, halfLengthA);

            sB = Math.Clamp(b * sA + e, -halfLengthB, halfLengthB);
            sA = Math.Clamp(b * sB - d, -halfLengthA, halfLengthA);

            onA = centerA + dirA * sA;
            onB = centerB + dirB * sB;
        }

        private static double Component(Vector3D v, int index)
            => index == 0 ? v.X : index == 1 ? v.Y : v.Z;

        private static double SignNonZero(double value) => value < 0 ? -1.0 : 1.0;

        private static bool TestAxis(
            Vector3D axis,
            Vector3D halfA, Vector3D halfB,
            Vector3D[] axesA, Vector3D[] axesB,
            Vector3D t,
            ref double minPenetration,
            ref Vector3D minAxis,
            ref int minAxisType,
            int axisType)
        {
            double projA = ProjectBox(halfA, axesA, axis);
            double projB = ProjectBox(halfB, axesB, axis);
            double distance = Math.Abs(Vector3D.Dot(t, axis));

            double penetration = projA + projB - distance;

            if (penetration < 0)
                return false;

            if (penetration < minPenetration)
            {
                minPenetration = penetration;
                minAxis = axis;
                minAxisType = axisType;
            }

            return true;
        }

        private static double ProjectBox(Vector3D halfExtents, Vector3D[] axes, Vector3D axis)
        {
            return halfExtents.X * Math.Abs(Vector3D.Dot(axes[0], axis)) +
                   halfExtents.Y * Math.Abs(Vector3D.Dot(axes[1], axis)) +
                   halfExtents.Z * Math.Abs(Vector3D.Dot(axes[2], axis));
        }

        #endregion

        #region Distance Queries

        /// <summary>
        /// Distance from a point to the surface of a body's collision shape.
        /// Returns 0 when the point is inside the shape.
        /// </summary>
        public static double DistanceToBody(Vector3D point, IPhysicsBody body)
        {
            if (body is RigidBody rb)
            {
                switch (rb.ShapeType)
                {
                    case CollisionShapeType.Sphere:
                        return Math.Max(0, (point - rb.Position).Magnitude - rb.Radius);

                    case CollisionShapeType.Box:
                        var local = rb.Transform.InverseTransformPoint(point);
                        var h = rb.HalfExtents;
                        double dx = Math.Max(0, Math.Abs(local.X) - h.X);
                        double dy = Math.Max(0, Math.Abs(local.Y) - h.Y);
                        double dz = Math.Max(0, Math.Abs(local.Z) - h.Z);
                        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
                }
            }

            return Math.Sqrt(body.BoundingBox.DistanceSquared(point));
        }

        #endregion

        #region Ray Casting

        /// <summary>
        /// Casts a ray against a rigid body.
        /// </summary>
        /// <param name="origin">Ray origin.</param>
        /// <param name="direction">Ray direction (normalized).</param>
        /// <param name="body">The body to test against.</param>
        /// <param name="maxDistance">Maximum ray distance.</param>
        /// <param name="hitDistance">Distance to hit point (if hit).</param>
        /// <param name="hitPoint">Hit point in world space (if hit).</param>
        /// <param name="hitNormal">Surface normal at hit point (if hit).</param>
        /// <returns>True if the ray hits the body.</returns>
        public static bool Raycast(
            Vector3D origin,
            Vector3D direction,
            RigidBody body,
            double maxDistance,
            out double hitDistance,
            out Vector3D hitPoint,
            out Vector3D hitNormal)
        {
            hitDistance = 0;
            hitPoint = Vector3D.Zero;
            hitNormal = Vector3D.Zero;

            switch (body.ShapeType)
            {
                case CollisionShapeType.Sphere:
                    return RaycastSphere(origin, direction, body.Position, body.Radius, maxDistance,
                        out hitDistance, out hitPoint, out hitNormal);

                case CollisionShapeType.Box:
                    return RaycastBox(origin, direction, body, maxDistance,
                        out hitDistance, out hitPoint, out hitNormal);

                default:
                    return false;
            }
        }

        private static bool RaycastSphere(
            Vector3D origin,
            Vector3D direction,
            Vector3D center,
            double radius,
            double maxDistance,
            out double hitDistance,
            out Vector3D hitPoint,
            out Vector3D hitNormal)
        {
            hitDistance = 0;
            hitPoint = Vector3D.Zero;
            hitNormal = Vector3D.Zero;

            var oc = origin - center;
            double a = Vector3D.Dot(direction, direction);
            double b = 2.0 * Vector3D.Dot(oc, direction);
            double c = Vector3D.Dot(oc, oc) - radius * radius;
            double discriminant = b * b - 4 * a * c;

            if (discriminant < 0)
                return false;

            double sqrtD = Math.Sqrt(discriminant);
            double t = (-b - sqrtD) / (2 * a);

            if (t < 0)
            {
                t = (-b + sqrtD) / (2 * a);
                if (t < 0)
                    return false;
            }

            if (t > maxDistance)
                return false;

            hitDistance = t;
            hitPoint = origin + direction * t;
            hitNormal = (hitPoint - center).Normalized;
            return true;
        }

        private static bool RaycastBox(
            Vector3D origin,
            Vector3D direction,
            RigidBody box,
            double maxDistance,
            out double hitDistance,
            out Vector3D hitPoint,
            out Vector3D hitNormal)
        {
            hitDistance = 0;
            hitPoint = Vector3D.Zero;
            hitNormal = Vector3D.Zero;

            // Transform ray to box local space
            var localOrigin = box.Transform.InverseTransformPoint(origin);
            var localDir = box.Transform.InverseTransformDirection(direction);

            var halfExtents = box.HalfExtents;
            var localAABB = new AABB(-halfExtents, halfExtents);

            if (!localAABB.RayIntersect(localOrigin, localDir, out double tMin, out double tMax))
                return false;

            if (tMin < 0)
                tMin = tMax;

            if (tMin < 0 || tMin > maxDistance)
                return false;

            hitDistance = tMin;
            var localHit = localOrigin + localDir * tMin;

            // Calculate local normal
            var localNormal = Vector3D.Zero;
            double bias = 0.0001;
            if (Math.Abs(localHit.X - halfExtents.X) < bias) localNormal = Vector3D.Right;
            else if (Math.Abs(localHit.X + halfExtents.X) < bias) localNormal = Vector3D.Left;
            else if (Math.Abs(localHit.Y - halfExtents.Y) < bias) localNormal = Vector3D.Up;
            else if (Math.Abs(localHit.Y + halfExtents.Y) < bias) localNormal = Vector3D.Down;
            else if (Math.Abs(localHit.Z - halfExtents.Z) < bias) localNormal = Vector3D.Forward;
            else if (Math.Abs(localHit.Z + halfExtents.Z) < bias) localNormal = Vector3D.Backward;

            hitPoint = box.Transform.TransformPoint(localHit);
            hitNormal = box.Transform.TransformDirection(localNormal);
            return true;
        }

        #endregion
    }
}
