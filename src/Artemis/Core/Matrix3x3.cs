using System;
using System.Runtime.CompilerServices;

namespace Artemis.Core
{
    /// <summary>
    /// Represents a 3x3 matrix for rotations, scaling, and inertia tensors.
    /// </summary>
    public struct Matrix3x3 : IEquatable<Matrix3x3>
    {
        public double M00, M01, M02;
        public double M10, M11, M12;
        public double M20, M21, M22;

        #region Static Matrices

        /// <summary>Returns the identity matrix.</summary>
        public static readonly Matrix3x3 Identity = new(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        );

        /// <summary>Returns a zero matrix.</summary>
        public static readonly Matrix3x3 Zero = new(
            0, 0, 0,
            0, 0, 0,
            0, 0, 0
        );

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new matrix with the specified components.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Matrix3x3(
            double m00, double m01, double m02,
            double m10, double m11, double m12,
            double m20, double m21, double m22)
        {
            M00 = m00; M01 = m01; M02 = m02;
            M10 = m10; M11 = m11; M12 = m12;
            M20 = m20; M21 = m21; M22 = m22;
        }

        /// <summary>
        /// Creates a diagonal matrix.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 Diagonal(double d0, double d1, double d2)
            => new(d0, 0, 0, 0, d1, 0, 0, 0, d2);

        /// <summary>
        /// Creates a diagonal matrix with the same value.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 Diagonal(double d)
            => new(d, 0, 0, 0, d, 0, 0, 0, d);

        /// <summary>
        /// Creates a matrix from column vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 FromColumns(Vector3D c0, Vector3D c1, Vector3D c2)
            => new(
                c0.X, c1.X, c2.X,
                c0.Y, c1.Y, c2.Y,
                c0.Z, c1.Z, c2.Z
            );

        /// <summary>
        /// Creates a matrix from row vectors.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 FromRows(Vector3D r0, Vector3D r1, Vector3D r2)
            => new(
                r0.X, r0.Y, r0.Z,
                r1.X, r1.Y, r1.Z,
                r2.X, r2.Y, r2.Z
            );

        /// <summary>
        /// Creates a rotation matrix from a quaternion.
        /// </summary>
        public static Matrix3x3 FromQuaternion(Quaternion q)
        {
            double x2 = q.X * 2, y2 = q.Y * 2, z2 = q.Z * 2;
            double xx2 = q.X * x2, yy2 = q.Y * y2, zz2 = q.Z * z2;
            double xy2 = q.X * y2, xz2 = q.X * z2, yz2 = q.Y * z2;
            double wx2 = q.W * x2, wy2 = q.W * y2, wz2 = q.W * z2;

            return new Matrix3x3(
                1 - yy2 - zz2, xy2 - wz2, xz2 + wy2,
                xy2 + wz2, 1 - xx2 - zz2, yz2 - wx2,
                xz2 - wy2, yz2 + wx2, 1 - xx2 - yy2
            );
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets the transpose of this matrix.
        /// </summary>
        public readonly Matrix3x3 Transposed
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => new(
                M00, M10, M20,
                M01, M11, M21,
                M02, M12, M22
            );
        }

        /// <summary>
        /// Gets the determinant of this matrix.
        /// </summary>
        public readonly double Determinant
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => M00 * (M11 * M22 - M12 * M21)
                 - M01 * (M10 * M22 - M12 * M20)
                 + M02 * (M10 * M21 - M11 * M20);
        }

        /// <summary>
        /// Gets the inverse of this matrix.
        /// </summary>
        public readonly Matrix3x3 Inverse
        {
            get
            {
                double det = Determinant;
                if (Math.Abs(det) < PhysicsConstants.Epsilon)
                    return Identity;

                double invDet = 1.0 / det;

                return new Matrix3x3(
                    (M11 * M22 - M12 * M21) * invDet,
                    (M02 * M21 - M01 * M22) * invDet,
                    (M01 * M12 - M02 * M11) * invDet,
                    (M12 * M20 - M10 * M22) * invDet,
                    (M00 * M22 - M02 * M20) * invDet,
                    (M02 * M10 - M00 * M12) * invDet,
                    (M10 * M21 - M11 * M20) * invDet,
                    (M01 * M20 - M00 * M21) * invDet,
                    (M00 * M11 - M01 * M10) * invDet
                );
            }
        }

        #endregion

        #region Operators

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 operator +(Matrix3x3 a, Matrix3x3 b)
            => new(
                a.M00 + b.M00, a.M01 + b.M01, a.M02 + b.M02,
                a.M10 + b.M10, a.M11 + b.M11, a.M12 + b.M12,
                a.M20 + b.M20, a.M21 + b.M21, a.M22 + b.M22
            );

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 operator -(Matrix3x3 a, Matrix3x3 b)
            => new(
                a.M00 - b.M00, a.M01 - b.M01, a.M02 - b.M02,
                a.M10 - b.M10, a.M11 - b.M11, a.M12 - b.M12,
                a.M20 - b.M20, a.M21 - b.M21, a.M22 - b.M22
            );

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 operator *(Matrix3x3 a, double s)
            => new(
                a.M00 * s, a.M01 * s, a.M02 * s,
                a.M10 * s, a.M11 * s, a.M12 * s,
                a.M20 * s, a.M21 * s, a.M22 * s
            );

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 operator *(double s, Matrix3x3 a)
            => a * s;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Matrix3x3 operator *(Matrix3x3 a, Matrix3x3 b)
            => new(
                a.M00 * b.M00 + a.M01 * b.M10 + a.M02 * b.M20,
                a.M00 * b.M01 + a.M01 * b.M11 + a.M02 * b.M21,
                a.M00 * b.M02 + a.M01 * b.M12 + a.M02 * b.M22,
                a.M10 * b.M00 + a.M11 * b.M10 + a.M12 * b.M20,
                a.M10 * b.M01 + a.M11 * b.M11 + a.M12 * b.M21,
                a.M10 * b.M02 + a.M11 * b.M12 + a.M12 * b.M22,
                a.M20 * b.M00 + a.M21 * b.M10 + a.M22 * b.M20,
                a.M20 * b.M01 + a.M21 * b.M11 + a.M22 * b.M21,
                a.M20 * b.M02 + a.M21 * b.M12 + a.M22 * b.M22
            );

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(Matrix3x3 m, Vector3D v)
            => new(
                m.M00 * v.X + m.M01 * v.Y + m.M02 * v.Z,
                m.M10 * v.X + m.M11 * v.Y + m.M12 * v.Z,
                m.M20 * v.X + m.M21 * v.Y + m.M22 * v.Z
            );

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Matrix3x3 a, Matrix3x3 b)
            => a.Equals(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Matrix3x3 a, Matrix3x3 b)
            => !a.Equals(b);

        #endregion

        #region Indexer

        /// <summary>
        /// Gets or sets the element at the specified row and column.
        /// </summary>
        public double this[int row, int col]
        {
            readonly get => (row, col) switch
            {
                (0, 0) => M00, (0, 1) => M01, (0, 2) => M02,
                (1, 0) => M10, (1, 1) => M11, (1, 2) => M12,
                (2, 0) => M20, (2, 1) => M21, (2, 2) => M22,
                _ => throw new IndexOutOfRangeException()
            };
            set
            {
                switch ((row, col))
                {
                    case (0, 0): M00 = value; break;
                    case (0, 1): M01 = value; break;
                    case (0, 2): M02 = value; break;
                    case (1, 0): M10 = value; break;
                    case (1, 1): M11 = value; break;
                    case (1, 2): M12 = value; break;
                    case (2, 0): M20 = value; break;
                    case (2, 1): M21 = value; break;
                    case (2, 2): M22 = value; break;
                    default: throw new IndexOutOfRangeException();
                }
            }
        }

        #endregion

        #region Equality

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool Equals(Matrix3x3 other)
            => Math.Abs(M00 - other.M00) < PhysicsConstants.Epsilon &&
               Math.Abs(M01 - other.M01) < PhysicsConstants.Epsilon &&
               Math.Abs(M02 - other.M02) < PhysicsConstants.Epsilon &&
               Math.Abs(M10 - other.M10) < PhysicsConstants.Epsilon &&
               Math.Abs(M11 - other.M11) < PhysicsConstants.Epsilon &&
               Math.Abs(M12 - other.M12) < PhysicsConstants.Epsilon &&
               Math.Abs(M20 - other.M20) < PhysicsConstants.Epsilon &&
               Math.Abs(M21 - other.M21) < PhysicsConstants.Epsilon &&
               Math.Abs(M22 - other.M22) < PhysicsConstants.Epsilon;

        public override readonly bool Equals(object? obj)
            => obj is Matrix3x3 other && Equals(other);

        public override readonly int GetHashCode()
            => HashCode.Combine(
                HashCode.Combine(M00, M01, M02),
                HashCode.Combine(M10, M11, M12),
                HashCode.Combine(M20, M21, M22)
            );

        public override readonly string ToString()
            => $"[({M00:F3}, {M01:F3}, {M02:F3}), ({M10:F3}, {M11:F3}, {M12:F3}), ({M20:F3}, {M21:F3}, {M22:F3})]";

        #endregion
    }
}
