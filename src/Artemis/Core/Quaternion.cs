using System;
using System.Runtime.CompilerServices;

namespace Artemis.Core
{
    /// <summary>
    /// Represents a quaternion used for 3D rotations.
    /// </summary>
    public struct Quaternion : IEquatable<Quaternion>
    {
        /// <summary>X component.</summary>
        public double X;

        /// <summary>Y component.</summary>
        public double Y;

        /// <summary>Z component.</summary>
        public double Z;

        /// <summary>W component (scalar part).</summary>
        public double W;

        #region Static Quaternions

        /// <summary>Returns the identity quaternion (no rotation).</summary>
        public static readonly Quaternion Identity = new(0, 0, 0, 1);

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new quaternion with the specified components.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Quaternion(double x, double y, double z, double w)
        {
            X = x;
            Y = y;
            Z = z;
            W = w;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets the magnitude of the quaternion.
        /// </summary>
        public readonly double Magnitude
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Math.Sqrt(X * X + Y * Y + Z * Z + W * W);
        }

        /// <summary>
        /// Gets the normalized version of this quaternion.
        /// </summary>
        public readonly Quaternion Normalized
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                double mag = Magnitude;
                if (mag > PhysicsConstants.Epsilon)
                    return new Quaternion(X / mag, Y / mag, Z / mag, W / mag);
                return Identity;
            }
        }


        /// <summary>
        /// Gets the conjugate of this quaternion.
        /// </summary>
        public readonly Quaternion Conjugate
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => new(-X, -Y, -Z, W);
        }

        /// <summary>
        /// Gets the inverse of this quaternion.
        /// </summary>
        public readonly Quaternion Inverse
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get
            {
                double sqrMag = X * X + Y * Y + Z * Z + W * W;
                if (sqrMag > PhysicsConstants.Epsilon)
                    return new Quaternion(-X / sqrMag, -Y / sqrMag, -Z / sqrMag, W / sqrMag);
                return Identity;
            }
        }

        #endregion

        #region Static Factory Methods

        /// <summary>
        /// Creates a quaternion from Euler angles (in radians).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion FromEuler(double pitch, double yaw, double roll)
        {
            double cy = Math.Cos(yaw * 0.5);
            double sy = Math.Sin(yaw * 0.5);
            double cp = Math.Cos(pitch * 0.5);
            double sp = Math.Sin(pitch * 0.5);
            double cr = Math.Cos(roll * 0.5);
            double sr = Math.Sin(roll * 0.5);

            return new Quaternion(
                sr * cp * cy - cr * sp * sy,
                cr * sp * cy + sr * cp * sy,
                cr * cp * sy - sr * sp * cy,
                cr * cp * cy + sr * sp * sy
            );
        }

        /// <summary>
        /// Creates a quaternion from Euler angles vector (in radians).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion FromEuler(Vector3D eulerAngles)
            => FromEuler(eulerAngles.X, eulerAngles.Y, eulerAngles.Z);

        /// <summary>
        /// Creates a quaternion representing a rotation around an axis.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion FromAxisAngle(Vector3D axis, double angle)
        {
            axis = axis.Normalized;
            double halfAngle = angle * 0.5;
            double sinHalf = Math.Sin(halfAngle);
            return new Quaternion(
                axis.X * sinHalf,
                axis.Y * sinHalf,
                axis.Z * sinHalf,
                Math.Cos(halfAngle)
            );
        }

        /// <summary>
        /// Creates a quaternion that rotates from one direction to another.
        /// </summary>
        public static Quaternion FromToRotation(Vector3D from, Vector3D to)
        {
            from = from.Normalized;
            to = to.Normalized;

            double dot = Vector3D.Dot(from, to);

            if (dot > 0.999999)
                return Identity;

            if (dot < -0.999999)
            {
                var axis = Vector3D.Cross(Vector3D.Right, from);
                if (axis.MagnitudeSquared < PhysicsConstants.Epsilon)
                    axis = Vector3D.Cross(Vector3D.Up, from);
                return FromAxisAngle(axis.Normalized, PhysicsConstants.Pi);
            }

            var cross = Vector3D.Cross(from, to);
            return new Quaternion(cross.X, cross.Y, cross.Z, 1 + dot).Normalized;
        }

        /// <summary>
        /// Creates a look rotation quaternion.
        /// </summary>
        public static Quaternion LookRotation(Vector3D forward, Vector3D up)
        {
            forward = forward.Normalized;
            var right = Vector3D.Cross(up, forward).Normalized;
            up = Vector3D.Cross(forward, right);

            double m00 = right.X, m01 = right.Y, m02 = right.Z;
            double m10 = up.X, m11 = up.Y, m12 = up.Z;
            double m20 = forward.X, m21 = forward.Y, m22 = forward.Z;

            double trace = m00 + m11 + m22;
            Quaternion q;

            if (trace > 0)
            {
                double s = 0.5 / Math.Sqrt(trace + 1.0);
                q = new Quaternion(
                    (m12 - m21) * s,
                    (m20 - m02) * s,
                    (m01 - m10) * s,
                    0.25 / s
                );
            }
            else if (m00 > m11 && m00 > m22)
            {
                double s = 2.0 * Math.Sqrt(1.0 + m00 - m11 - m22);
                q = new Quaternion(
                    0.25 * s,
                    (m01 + m10) / s,
                    (m02 + m20) / s,
                    (m12 - m21) / s
                );
            }
            else if (m11 > m22)
            {
                double s = 2.0 * Math.Sqrt(1.0 + m11 - m00 - m22);
                q = new Quaternion(
                    (m01 + m10) / s,
                    0.25 * s,
                    (m12 + m21) / s,
                    (m20 - m02) / s
                );
            }
            else
            {
                double s = 2.0 * Math.Sqrt(1.0 + m22 - m00 - m11);
                q = new Quaternion(
                    (m02 + m20) / s,
                    (m12 + m21) / s,
                    0.25 * s,
                    (m01 - m10) / s
                );
            }

            return q.Normalized;
        }

        #endregion

        #region Operators

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion operator *(Quaternion a, Quaternion b)
            => new(
                a.W * b.X + a.X * b.W + a.Y * b.Z - a.Z * b.Y,
                a.W * b.Y - a.X * b.Z + a.Y * b.W + a.Z * b.X,
                a.W * b.Z + a.X * b.Y - a.Y * b.X + a.Z * b.W,
                a.W * b.W - a.X * b.X - a.Y * b.Y - a.Z * b.Z
            );

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion operator +(Quaternion a, Quaternion b)
            => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z, a.W + b.W);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion operator -(Quaternion a, Quaternion b)
            => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z, a.W - b.W);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion operator *(Quaternion q, double scalar)
            => new(q.X * scalar, q.Y * scalar, q.Z * scalar, q.W * scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion operator *(double scalar, Quaternion q)
            => q * scalar;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Quaternion operator /(Quaternion q, double scalar)
            => new(q.X / scalar, q.Y / scalar, q.Z / scalar, q.W / scalar);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(Quaternion q, Vector3D v)
        {
            double x2 = q.X * 2, y2 = q.Y * 2, z2 = q.Z * 2;
            double xx2 = q.X * x2, yy2 = q.Y * y2, zz2 = q.Z * z2;
            double xy2 = q.X * y2, xz2 = q.X * z2, yz2 = q.Y * z2;
            double wx2 = q.W * x2, wy2 = q.W * y2, wz2 = q.W * z2;

            return new Vector3D(
                (1 - yy2 - zz2) * v.X + (xy2 - wz2) * v.Y + (xz2 + wy2) * v.Z,
                (xy2 + wz2) * v.X + (1 - xx2 - zz2) * v.Y + (yz2 - wx2) * v.Z,
                (xz2 - wy2) * v.X + (yz2 + wx2) * v.Y + (1 - xx2 - yy2) * v.Z
            );
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Quaternion a, Quaternion b)
            => a.Equals(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Quaternion a, Quaternion b)
            => !a.Equals(b);

        #endregion

        #region Static Methods

        /// <summary>
        /// Spherically interpolates between two quaternions.
        /// </summary>
        public static Quaternion Slerp(Quaternion a, Quaternion b, double t)
        {
            t = Math.Clamp(t, 0.0, 1.0);

            double dot = a.X * b.X + a.Y * b.Y + a.Z * b.Z + a.W * b.W;

            if (dot < 0)
            {
                b = new Quaternion(-b.X, -b.Y, -b.Z, -b.W);
                dot = -dot;
            }

            if (dot > 0.9995)
            {
                return new Quaternion(
                    a.X + t * (b.X - a.X),
                    a.Y + t * (b.Y - a.Y),
                    a.Z + t * (b.Z - a.Z),
                    a.W + t * (b.W - a.W)
                ).Normalized;
            }

            double theta0 = Math.Acos(dot);
            double theta = theta0 * t;
            double sinTheta = Math.Sin(theta);
            double sinTheta0 = Math.Sin(theta0);

            double s0 = Math.Cos(theta) - dot * sinTheta / sinTheta0;
            double s1 = sinTheta / sinTheta0;

            return new Quaternion(
                s0 * a.X + s1 * b.X,
                s0 * a.Y + s1 * b.Y,
                s0 * a.Z + s1 * b.Z,
                s0 * a.W + s1 * b.W
            );
        }

        /// <summary>
        /// Returns the angle in radians between two quaternions.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Angle(Quaternion a, Quaternion b)
        {
            double dot = Math.Abs(a.X * b.X + a.Y * b.Y + a.Z * b.Z + a.W * b.W);
            return 2.0 * Math.Acos(Math.Min(dot, 1.0));
        }

        #endregion

        #region Instance Methods

        private static double CopySign(double value, double sign)
        {
#if NETSTANDARD2_1
            long valueBits = BitConverter.DoubleToInt64Bits(value);
            long signBits = BitConverter.DoubleToInt64Bits(sign);
            valueBits = (valueBits & 0x7FFFFFFFFFFFFFFF) | (signBits & unchecked((long)0x8000000000000000));
            return BitConverter.Int64BitsToDouble(valueBits);
#else
            return Math.CopySign(value, sign);
#endif
        }

        /// <summary>
        /// Converts this quaternion to Euler angles (in radians).
        /// </summary>
        public readonly Vector3D ToEuler()
        {
            double sinrCosp = 2 * (W * X + Y * Z);
            double cosrCosp = 1 - 2 * (X * X + Y * Y);
            double roll = Math.Atan2(sinrCosp, cosrCosp);

            double sinp = 2 * (W * Y - Z * X);
            double pitch = Math.Abs(sinp) >= 1
                ? CopySign(PhysicsConstants.Pi / 2, sinp)
                : Math.Asin(sinp);

            double sinyCosp = 2 * (W * Z + X * Y);
            double cosyCosp = 1 - 2 * (Y * Y + Z * Z);
            double yaw = Math.Atan2(sinyCosp, cosyCosp);

            return new Vector3D(pitch, yaw, roll);
        }

        /// <summary>
        /// Normalizes this quaternion in place.
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
                W /= mag;
            }
        }

        #endregion

        #region Equality

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Equals(Quaternion other)
            => Math.Abs(X - other.X) < PhysicsConstants.Epsilon &&
               Math.Abs(Y - other.Y) < PhysicsConstants.Epsilon &&
               Math.Abs(Z - other.Z) < PhysicsConstants.Epsilon &&
               Math.Abs(W - other.W) < PhysicsConstants.Epsilon;

        public override readonly bool Equals(object? obj)
            => obj is Quaternion other && Equals(other);

        public override readonly int GetHashCode()
            => HashCode.Combine(X, Y, Z, W);

        public override readonly string ToString()
            => $"({X:F3}, {Y:F3}, {Z:F3}, {W:F3})";

        #endregion
    }
}
