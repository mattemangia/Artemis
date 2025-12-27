using System;
using System.Runtime.CompilerServices;

namespace Artemis.Core
{
    /// <summary>
    /// Represents a 3D vector with double precision.
    /// Optimized for physics calculations and compatible with Unity/MonoGame conversions.
    /// </summary>
    public struct Vector3D : IEquatable<Vector3D>
    {
        /// <summary>X component of the vector.</summary>
        public double X;

        /// <summary>Y component of the vector.</summary>
        public double Y;

        /// <summary>Z component of the vector.</summary>
        public double Z;

        #region Static Vectors

        /// <summary>Returns a vector with all components set to zero.</summary>
        public static readonly Vector3D Zero = new(0, 0, 0);

        /// <summary>Returns a vector with all components set to one.</summary>
        public static readonly Vector3D One = new(1, 1, 1);

        /// <summary>Returns the up vector (0, 1, 0).</summary>
        public static readonly Vector3D Up = new(0, 1, 0);

        /// <summary>Returns the down vector (0, -1, 0).</summary>
        public static readonly Vector3D Down = new(0, -1, 0);

        /// <summary>Returns the right vector (1, 0, 0).</summary>
        public static readonly Vector3D Right = new(1, 0, 0);

        /// <summary>Returns the left vector (-1, 0, 0).</summary>
        public static readonly Vector3D Left = new(-1, 0, 0);

        /// <summary>Returns the forward vector (0, 0, 1).</summary>
        public static readonly Vector3D Forward = new(0, 0, 1);

        /// <summary>Returns the backward vector (0, 0, -1).</summary>
        public static readonly Vector3D Backward = new(0, 0, -1);

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new Vector3D with the specified components.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        /// <summary>
        /// Creates a new Vector3D with all components set to the same value.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector3D(double value)
        {
            X = Y = Z = value;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets the magnitude (length) of the vector.
        /// </summary>
        public readonly double Magnitude
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Math.Sqrt(X * X + Y * Y + Z * Z);
        }

        /// <summary>
        /// Gets the length (alias for Magnitude).
        /// </summary>
        public readonly double Length
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Magnitude;
        }

        /// <summary>
        /// Gets the squared magnitude of the vector (faster than Magnitude).
        /// </summary>
        public readonly double MagnitudeSquared
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => X * X + Y * Y + Z * Z;
        }

        /// <summary>
        /// Gets the squared length (alias for MagnitudeSquared).
        /// </summary>
        public readonly double LengthSquared
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => MagnitudeSquared;
        }

        /// <summary>
        /// Gets the normalized version of this vector.
        /// </summary>
        public readonly Vector3D Normalized
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                double mag = Magnitude;
                if (mag > PhysicsConstants.Epsilon)
                    return this / mag;
                return Zero;
            }
        }


        #endregion

        #region Operators

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator +(Vector3D a, Vector3D b)
            => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator -(Vector3D a, Vector3D b)
            => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(Vector3D v, double scalar)
            => new(v.X * scalar, v.Y * scalar, v.Z * scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(double scalar, Vector3D v)
            => new(v.X * scalar, v.Y * scalar, v.Z * scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator /(Vector3D v, double scalar)
            => new(v.X / scalar, v.Y / scalar, v.Z / scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator System.Numerics.Vector3(Vector3D v)
            => new((float)v.X, (float)v.Y, (float)v.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Vector3D(System.Numerics.Vector3 v)
            => new(v.X, v.Y, v.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator -(Vector3D v)
            => new(-v.X, -v.Y, -v.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Vector3D a, Vector3D b)
            => a.Equals(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Vector3D a, Vector3D b)
            => !a.Equals(b);

        #endregion

        #region Static Methods

        /// <summary>
        /// Calculates the dot product of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(Vector3D a, Vector3D b)
            => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

        /// <summary>
        /// Calculates the cross product of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Cross(Vector3D a, Vector3D b)
            => new(
                a.Y * b.Z - a.Z * b.Y,
                a.Z * b.X - a.X * b.Z,
                a.X * b.Y - a.Y * b.X
            );

        /// <summary>
        /// Calculates the distance between two points.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Distance(Vector3D a, Vector3D b)
            => (b - a).Magnitude;

        /// <summary>
        /// Calculates the squared distance between two points (faster than Distance).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double DistanceSquared(Vector3D a, Vector3D b)
            => (b - a).MagnitudeSquared;

        /// <summary>
        /// Linearly interpolates between two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Lerp(Vector3D a, Vector3D b, double t)
        {
            t = Math.Clamp(t, 0.0, 1.0);
            return new Vector3D(
                a.X + (b.X - a.X) * t,
                a.Y + (b.Y - a.Y) * t,
                a.Z + (b.Z - a.Z) * t
            );
        }

        /// <summary>
        /// Returns the component-wise minimum of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Min(Vector3D a, Vector3D b)
            => new(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y), Math.Min(a.Z, b.Z));

        /// <summary>
        /// Returns the component-wise maximum of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Max(Vector3D a, Vector3D b)
            => new(Math.Max(a.X, b.X), Math.Max(a.Y, b.Y), Math.Max(a.Z, b.Z));

        /// <summary>
        /// Reflects a vector off a surface with the specified normal.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Reflect(Vector3D vector, Vector3D normal)
            => vector - 2.0 * Dot(vector, normal) * normal;

        /// <summary>
        /// Projects a vector onto another vector.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Project(Vector3D vector, Vector3D onto)
        {
            double dot = Dot(onto, onto);
            if (dot < PhysicsConstants.Epsilon)
                return Zero;
            return onto * (Dot(vector, onto) / dot);
        }

        /// <summary>
        /// Returns the angle in radians between two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Angle(Vector3D from, Vector3D to)
        {
            double dot = Dot(from.Normalized, to.Normalized);
            return Math.Acos(Math.Clamp(dot, -1.0, 1.0));
        }

        #endregion

        #region Instance Methods

        /// <summary>
        /// Normalizes this vector in place.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Normalize()
        {
            double mag = Magnitude;
            if (mag > PhysicsConstants.Epsilon)
            {
                X /= mag;
                Y /= mag;
                Z /= mag;
            }
        }

        /// <summary>
        /// Clamps the magnitude of the vector.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D ClampMagnitude(double maxMagnitude)
        {
            double sqrMag = MagnitudeSquared;
            if (sqrMag > maxMagnitude * maxMagnitude)
            {
                double mag = Math.Sqrt(sqrMag);
                return this * (maxMagnitude / mag);
            }
            return this;
        }

        #endregion

        #region Equality & Comparison

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Equals(Vector3D other)
            => Math.Abs(X - other.X) < PhysicsConstants.Epsilon &&
               Math.Abs(Y - other.Y) < PhysicsConstants.Epsilon &&
               Math.Abs(Z - other.Z) < PhysicsConstants.Epsilon;

        public override readonly bool Equals(object? obj)
            => obj is Vector3D other && Equals(other);

        public override readonly int GetHashCode()
            => HashCode.Combine(X, Y, Z);

        public override readonly string ToString()
            => $"({X:F3}, {Y:F3}, {Z:F3})";

        #endregion
    }
}
