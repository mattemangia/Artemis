using System;
using System.Runtime.CompilerServices;

namespace Artemis.Core
{
    /// <summary>
    /// Represents an Axis-Aligned Bounding Box for fast collision detection.
    /// </summary>
    public struct AABB : IEquatable<AABB>
    {
        /// <summary>Minimum corner of the bounding box.</summary>
        public Vector3D Min;

        /// <summary>Maximum corner of the bounding box.</summary>
        public Vector3D Max;

        #region Constructors

        /// <summary>
        /// Creates a new AABB with the specified min and max corners.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public AABB(Vector3D min, Vector3D max)
        {
            Min = min;
            Max = max;
        }

        /// <summary>
        /// Creates a new AABB centered at a position with the specified size.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static AABB FromCenterSize(Vector3D center, Vector3D size)
        {
            var halfSize = size * 0.5;
            return new AABB(center - halfSize, center + halfSize);
        }

        /// <summary>
        /// Creates an AABB that contains all the specified points.
        /// </summary>
        public static AABB FromPoints(params Vector3D[] points)
        {
            if (points.Length == 0)
                return new AABB(Vector3D.Zero, Vector3D.Zero);

            var min = points[0];
            var max = points[0];

            for (int i = 1; i < points.Length; i++)
            {
                min = Vector3D.Min(min, points[i]);
                max = Vector3D.Max(max, points[i]);
            }

            return new AABB(min, max);
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets the center of the bounding box.
        /// </summary>
        public readonly Vector3D Center
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => (Min + Max) * 0.5;
        }

        /// <summary>
        /// Gets the size of the bounding box.
        /// </summary>
        public readonly Vector3D Size
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Max - Min;
        }

        /// <summary>
        /// Gets the half-size (extents) of the bounding box.
        /// </summary>
        public readonly Vector3D Extents
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => (Max - Min) * 0.5;
        }

        /// <summary>
        /// Gets the volume of the bounding box.
        /// </summary>
        public readonly double Volume
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                var size = Size;
                return size.X * size.Y * size.Z;
            }
        }

        /// <summary>
        /// Gets the surface area of the bounding box.
        /// </summary>
        public readonly double SurfaceArea
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                var size = Size;
                return 2.0 * (size.X * size.Y + size.Y * size.Z + size.Z * size.X);
            }
        }

        /// <summary>
        /// Gets the corners of the bounding box.
        /// </summary>
        public readonly Vector3D[] Corners
        {
            get
            {
                return new[]
                {
                    new Vector3D(Min.X, Min.Y, Min.Z),
                    new Vector3D(Max.X, Min.Y, Min.Z),
                    new Vector3D(Min.X, Max.Y, Min.Z),
                    new Vector3D(Max.X, Max.Y, Min.Z),
                    new Vector3D(Min.X, Min.Y, Max.Z),
                    new Vector3D(Max.X, Min.Y, Max.Z),
                    new Vector3D(Min.X, Max.Y, Max.Z),
                    new Vector3D(Max.X, Max.Y, Max.Z)
                };
            }
        }

        #endregion

        #region Methods

        /// <summary>
        /// Tests if this AABB contains a point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Contains(Vector3D point)
            => point.X >= Min.X && point.X <= Max.X &&
               point.Y >= Min.Y && point.Y <= Max.Y &&
               point.Z >= Min.Z && point.Z <= Max.Z;

        /// <summary>
        /// Tests if this AABB fully contains another AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Contains(AABB other)
            => Min.X <= other.Min.X && Max.X >= other.Max.X &&
               Min.Y <= other.Min.Y && Max.Y >= other.Max.Y &&
               Min.Z <= other.Min.Z && Max.Z >= other.Max.Z;

        /// <summary>
        /// Tests if this AABB intersects with another AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Intersects(AABB other)
            => Min.X <= other.Max.X && Max.X >= other.Min.X &&
               Min.Y <= other.Max.Y && Max.Y >= other.Min.Y &&
               Min.Z <= other.Max.Z && Max.Z >= other.Min.Z;

        /// <summary>
        /// Returns the union of this AABB and another AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly AABB Union(AABB other)
            => new(Vector3D.Min(Min, other.Min), Vector3D.Max(Max, other.Max));

        /// <summary>
        /// Returns the intersection of this AABB and another AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly AABB Intersection(AABB other)
            => new(Vector3D.Max(Min, other.Min), Vector3D.Min(Max, other.Max));

        /// <summary>
        /// Expands the AABB to include a point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Encapsulate(Vector3D point)
        {
            Min = Vector3D.Min(Min, point);
            Max = Vector3D.Max(Max, point);
        }

        /// <summary>
        /// Expands the AABB to include another AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Encapsulate(AABB other)
        {
            Min = Vector3D.Min(Min, other.Min);
            Max = Vector3D.Max(Max, other.Max);
        }

        /// <summary>
        /// Expands the AABB by the specified amount in all directions.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly AABB Expanded(double amount)
            => new(Min - new Vector3D(amount), Max + new Vector3D(amount));

        /// <summary>
        /// Gets the closest point on or inside the AABB to the specified point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D ClosestPoint(Vector3D point)
            => new(
                Math.Clamp(point.X, Min.X, Max.X),
                Math.Clamp(point.Y, Min.Y, Max.Y),
                Math.Clamp(point.Z, Min.Z, Max.Z)
            );

        /// <summary>
        /// Gets the squared distance from a point to the AABB.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly double DistanceSquared(Vector3D point)
        {
            var closest = ClosestPoint(point);
            return Vector3D.DistanceSquared(point, closest);
        }

        /// <summary>
        /// Ray-AABB intersection test.
        /// </summary>
        /// <param name="origin">Ray origin.</param>
        /// <param name="direction">Ray direction (should be normalized).</param>
        /// <param name="tMin">Minimum distance along ray (output).</param>
        /// <param name="tMax">Maximum distance along ray (output).</param>
        /// <returns>True if the ray intersects the AABB.</returns>
        public readonly bool RayIntersect(Vector3D origin, Vector3D direction, out double tMin, out double tMax)
        {
            tMin = double.MinValue;
            tMax = double.MaxValue;

            for (int i = 0; i < 3; i++)
            {
                double minVal = i == 0 ? Min.X : (i == 1 ? Min.Y : Min.Z);
                double maxVal = i == 0 ? Max.X : (i == 1 ? Max.Y : Max.Z);
                double originVal = i == 0 ? origin.X : (i == 1 ? origin.Y : origin.Z);
                double dirVal = i == 0 ? direction.X : (i == 1 ? direction.Y : direction.Z);

                if (Math.Abs(dirVal) < PhysicsConstants.Epsilon)
                {
                    if (originVal < minVal || originVal > maxVal)
                        return false;
                }
                else
                {
                    double invD = 1.0 / dirVal;
                    double t1 = (minVal - originVal) * invD;
                    double t2 = (maxVal - originVal) * invD;

                    if (t1 > t2)
                        (t1, t2) = (t2, t1);

                    tMin = Math.Max(tMin, t1);
                    tMax = Math.Min(tMax, t2);

                    if (tMin > tMax)
                        return false;
                }
            }

            return true;
        }

        #endregion

        #region Equality

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Equals(AABB other)
            => Min == other.Min && Max == other.Max;

        public override readonly bool Equals(object? obj)
            => obj is AABB other && Equals(other);

        public override readonly int GetHashCode()
            => HashCode.Combine(Min, Max);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(AABB a, AABB b)
            => a.Equals(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(AABB a, AABB b)
            => !a.Equals(b);

        public override readonly string ToString()
            => $"AABB(Min: {Min}, Max: {Max})";

        #endregion
    }
}
