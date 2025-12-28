using System;
using System.Runtime.CompilerServices;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Represents a 2D vector with double precision.
    /// Optimized for 2D physics calculations and compatible with Unity/MonoGame conversions.
    /// </summary>
    public struct Vector2D : IEquatable<Vector2D>
    {
        /// <summary>X component of the vector.</summary>
        public double X;

        /// <summary>Y component of the vector.</summary>
        public double Y;

        #region Static Vectors

        /// <summary>Returns a vector with all components set to zero.</summary>
        public static readonly Vector2D Zero = new(0, 0);

        /// <summary>Returns a vector with all components set to one.</summary>
        public static readonly Vector2D One = new(1, 1);

        /// <summary>Returns the up vector (0, 1).</summary>
        public static readonly Vector2D Up = new(0, 1);

        /// <summary>Returns the down vector (0, -1).</summary>
        public static readonly Vector2D Down = new(0, -1);

        /// <summary>Returns the right vector (1, 0).</summary>
        public static readonly Vector2D Right = new(1, 0);

        /// <summary>Returns the left vector (-1, 0).</summary>
        public static readonly Vector2D Left = new(-1, 0);

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new Vector2D with the specified components.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        /// Creates a new Vector2D with both components set to the same value.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D(double value)
        {
            X = Y = value;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets the magnitude (length) of the vector.
        /// </summary>
        public readonly double Magnitude
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Math.Sqrt(X * X + Y * Y);
        }

        /// <summary>
        /// Gets the magnitude (length) of the vector. Alias for Magnitude.
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
            get => X * X + Y * Y;
        }

        /// <summary>
        /// Gets the squared magnitude of the vector. Alias for MagnitudeSquared.
        /// </summary>
        public readonly double LengthSquared
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => MagnitudeSquared;
        }

        /// <summary>
        /// Gets the normalized version of this vector.
        /// </summary>
        public readonly Vector2D Normalized
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

        /// <summary>
        /// Gets the perpendicular vector (rotated 90 degrees counter-clockwise).
        /// </summary>
        public readonly Vector2D Perpendicular
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => new(-Y, X);
        }

        /// <summary>
        /// Gets the angle of this vector in radians from the positive X axis.
        /// </summary>
        public readonly double Angle
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Math.Atan2(Y, X);
        }

        #endregion

        #region Operators

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D operator +(Vector2D a, Vector2D b)
            => new(a.X + b.X, a.Y + b.Y);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D operator -(Vector2D a, Vector2D b)
            => new(a.X - b.X, a.Y - b.Y);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D operator *(Vector2D v, double scalar)
            => new(v.X * scalar, v.Y * scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D operator *(double scalar, Vector2D v)
            => new(v.X * scalar, v.Y * scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D operator /(Vector2D v, double scalar)
            => new(v.X / scalar, v.Y / scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D operator -(Vector2D v)
            => new(-v.X, -v.Y);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Vector2D a, Vector2D b)
            => a.Equals(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Vector2D a, Vector2D b)
            => !a.Equals(b);

        #endregion

        #region Static Methods

        /// <summary>
        /// Calculates the dot product of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(Vector2D a, Vector2D b)
            => a.X * b.X + a.Y * b.Y;

        /// <summary>
        /// Calculates the cross product (returns scalar in 2D - the Z component).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Cross(Vector2D a, Vector2D b)
            => a.X * b.Y - a.Y * b.X;

        /// <summary>
        /// Creates a vector perpendicular to the input, scaled by the scalar.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Cross(double s, Vector2D v)
            => new(-s * v.Y, s * v.X);

        /// <summary>
        /// Creates a vector perpendicular to the input, scaled by the scalar.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Cross(Vector2D v, double s)
            => new(s * v.Y, -s * v.X);

        /// <summary>
        /// Calculates the distance between two points.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Distance(Vector2D a, Vector2D b)
            => (b - a).Magnitude;

        /// <summary>
        /// Calculates the squared distance between two points (faster than Distance).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double DistanceSquared(Vector2D a, Vector2D b)
            => (b - a).MagnitudeSquared;

        /// <summary>
        /// Linearly interpolates between two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Lerp(Vector2D a, Vector2D b, double t)
        {
            t = Math.Clamp(t, 0.0, 1.0);
            return new Vector2D(
                a.X + (b.X - a.X) * t,
                a.Y + (b.Y - a.Y) * t
            );
        }

        /// <summary>
        /// Returns the component-wise minimum of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Min(Vector2D a, Vector2D b)
            => new(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y));

        /// <summary>
        /// Returns the component-wise maximum of two vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Max(Vector2D a, Vector2D b)
            => new(Math.Max(a.X, b.X), Math.Max(a.Y, b.Y));

        /// <summary>
        /// Reflects a vector off a surface with the specified normal.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Reflect(Vector2D vector, Vector2D normal)
            => vector - 2.0 * Dot(vector, normal) * normal;

        /// <summary>
        /// Projects a vector onto another vector.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Project(Vector2D vector, Vector2D onto)
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
        public static double AngleBetween(Vector2D from, Vector2D to)
        {
            double dot = Dot(from.Normalized, to.Normalized);
            return Math.Acos(Math.Clamp(dot, -1.0, 1.0));
        }

        /// <summary>
        /// Returns the signed angle in radians from one vector to another.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double SignedAngle(Vector2D from, Vector2D to)
        {
            double angle = AngleBetween(from, to);
            double cross = Cross(from, to);
            return cross < 0 ? -angle : angle;
        }

        /// <summary>
        /// Rotates a vector by the given angle in radians.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D Rotate(Vector2D v, double angleRadians)
        {
            double cos = Math.Cos(angleRadians);
            double sin = Math.Sin(angleRadians);
            return new Vector2D(
                v.X * cos - v.Y * sin,
                v.X * sin + v.Y * cos
            );
        }

        /// <summary>
        /// Creates a vector from an angle in radians.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D FromAngle(double angleRadians)
            => new(Math.Cos(angleRadians), Math.Sin(angleRadians));

        /// <summary>
        /// Creates a vector from an angle in radians with a specified length.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector2D FromAngle(double angleRadians, double length)
            => new(Math.Cos(angleRadians) * length, Math.Sin(angleRadians) * length);

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
            }
        }

        /// <summary>
        /// Clamps the magnitude of the vector.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector2D ClampMagnitude(double maxMagnitude)
        {
            double sqrMag = MagnitudeSquared;
            if (sqrMag > maxMagnitude * maxMagnitude)
            {
                double mag = Math.Sqrt(sqrMag);
                return this * (maxMagnitude / mag);
            }
            return this;
        }

        /// <summary>
        /// Converts to 3D vector with Z = 0.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D ToVector3D()
            => new(X, Y, 0);

        /// <summary>
        /// Converts to 3D vector with specified Z.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D ToVector3D(double z)
            => new(X, Y, z);

        #endregion

        #region Equality & Comparison

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Equals(Vector2D other)
            => Math.Abs(X - other.X) < PhysicsConstants.Epsilon &&
               Math.Abs(Y - other.Y) < PhysicsConstants.Epsilon;

        public override readonly bool Equals(object? obj)
            => obj is Vector2D other && Equals(other);

        public override readonly int GetHashCode()
            => HashCode.Combine(X, Y);

        public override readonly string ToString()
            => $"({X:F3}, {Y:F3})";

        #endregion
    }
}
