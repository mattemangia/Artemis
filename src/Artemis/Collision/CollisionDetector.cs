using System;
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
            bool inside = localCenter == closestLocal;

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
                        localNormal = new Vector3D(Math.Sign(localCenter.X), 0, 0);
                        closestLocal = new Vector3D(
                            Math.Sign(localCenter.X) * box.HalfExtents.X,
                            localCenter.Y,
                            localCenter.Z);
                        break;
                    case 1:
                        localNormal = new Vector3D(0, Math.Sign(localCenter.Y), 0);
                        closestLocal = new Vector3D(
                            localCenter.X,
                            Math.Sign(localCenter.Y) * box.HalfExtents.Y,
                            localCenter.Z);
                        break;
                    case 2:
                        localNormal = new Vector3D(0, 0, Math.Sign(localCenter.Z));
                        closestLocal = new Vector3D(
                            localCenter.X,
                            localCenter.Y,
                            Math.Sign(localCenter.Z) * box.HalfExtents.Z);
                        break;
                }

                normal = box.Transform.TransformDirection(localNormal);
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
                normal = box.Transform.TransformDirection(localNormal);
                penetration = sphere.Radius - distance;
                contactPoint = box.Transform.TransformPoint(closestLocal);
            }

            return CollisionInfo.Create(sphere, box, normal, penetration, contactPoint);
        }

        #endregion

        #region Box-Box (SAT)

        /// <summary>
        /// Detects collision between two boxes using the Separating Axis Theorem.
        /// </summary>
        public static CollisionInfo BoxBox(RigidBody boxA, RigidBody boxB)
        {
            // Get box orientations
            var rotA = Matrix3x3.FromQuaternion(boxA.Rotation);
            var rotB = Matrix3x3.FromQuaternion(boxB.Rotation);

            // Get axes for both boxes
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

            double minPenetration = double.MaxValue;
            Vector3D minAxis = Vector3D.Zero;
            int minAxisType = -1;

            // Test face axes of A
            for (int i = 0; i < 3; i++)
            {
                if (!TestAxis(axesA[i], halfA, halfB, axesA, axesB, t,
                    ref minPenetration, ref minAxis, ref minAxisType, i))
                    return CollisionInfo.NoCollision;
            }

            // Test face axes of B
            for (int i = 0; i < 3; i++)
            {
                if (!TestAxis(axesB[i], halfA, halfB, axesA, axesB, t,
                    ref minPenetration, ref minAxis, ref minAxisType, i + 3))
                    return CollisionInfo.NoCollision;
            }

            // Test edge axes (cross products)
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    var axis = Vector3D.Cross(axesA[i], axesB[j]);
                    if (axis.MagnitudeSquared < PhysicsConstants.Epsilon)
                        continue;

                    axis = axis.Normalized;
                    if (!TestAxis(axis, halfA, halfB, axesA, axesB, t,
                        ref minPenetration, ref minAxis, ref minAxisType, 6 + i * 3 + j))
                        return CollisionInfo.NoCollision;
                }
            }

            // Ensure normal points from A to B
            if (Vector3D.Dot(minAxis, t) < 0)
                minAxis = -minAxis;

            // Calculate contact point (simplified - center of overlap region)
            var contactPoint = (boxA.Position + boxB.Position) * 0.5;

            return CollisionInfo.Create(boxA, boxB, minAxis, minPenetration, contactPoint);
        }

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
