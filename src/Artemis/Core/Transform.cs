using System;
using System.Runtime.CompilerServices;

namespace Artemis.Core
{
    /// <summary>
    /// Represents a 3D transformation with position, rotation, and scale.
    /// </summary>
    public struct Transform : IEquatable<Transform>
    {
        /// <summary>Position in world space.</summary>
        public Vector3D Position;

        /// <summary>Rotation as a quaternion.</summary>
        public Quaternion Rotation;

        /// <summary>Scale factor.</summary>
        public Vector3D Scale;

        #region Static Transforms

        /// <summary>Returns the identity transform.</summary>
        public static readonly Transform Identity = new(Vector3D.Zero, Quaternion.Identity, Vector3D.One);

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new transform with the specified position, rotation, and scale.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Transform(Vector3D position, Quaternion rotation, Vector3D scale)
        {
            Position = position;
            Rotation = rotation;
            Scale = scale;
        }

        /// <summary>
        /// Creates a new transform with the specified position.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Transform(Vector3D position)
        {
            Position = position;
            Rotation = Quaternion.Identity;
            Scale = Vector3D.One;
        }

        /// <summary>
        /// Creates a new transform with the specified position and rotation.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Transform(Vector3D position, Quaternion rotation)
        {
            Position = position;
            Rotation = rotation;
            Scale = Vector3D.One;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets the forward direction of this transform.
        /// </summary>
        public readonly Vector3D Forward
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Rotation * Vector3D.Forward;
        }

        /// <summary>
        /// Gets the right direction of this transform.
        /// </summary>
        public readonly Vector3D Right
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Rotation * Vector3D.Right;
        }

        /// <summary>
        /// Gets the up direction of this transform.
        /// </summary>
        public readonly Vector3D Up
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Rotation * Vector3D.Up;
        }

        /// <summary>
        /// Gets the rotation matrix for this transform.
        /// </summary>
        public readonly Matrix3x3 RotationMatrix
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Matrix3x3.FromQuaternion(Rotation);
        }

        #endregion

        #region Methods

        /// <summary>
        /// Transforms a point from local to world space.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D TransformPoint(Vector3D localPoint)
        {
            var scaled = new Vector3D(
                localPoint.X * Scale.X,
                localPoint.Y * Scale.Y,
                localPoint.Z * Scale.Z
            );
            return Position + Rotation * scaled;
        }

        /// <summary>
        /// Transforms a direction from local to world space (ignores position and scale).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D TransformDirection(Vector3D localDirection)
            => Rotation * localDirection;

        /// <summary>
        /// Transforms a vector from local to world space (applies rotation and scale, ignores position).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D TransformVector(Vector3D localVector)
        {
            var scaled = new Vector3D(
                localVector.X * Scale.X,
                localVector.Y * Scale.Y,
                localVector.Z * Scale.Z
            );
            return Rotation * scaled;
        }

        /// <summary>
        /// Transforms a point from world to local space.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D InverseTransformPoint(Vector3D worldPoint)
        {
            var relative = worldPoint - Position;
            var rotated = Rotation.Inverse * relative;
            return new Vector3D(
                rotated.X / Scale.X,
                rotated.Y / Scale.Y,
                rotated.Z / Scale.Z
            );
        }

        /// <summary>
        /// Transforms a direction from world to local space.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly Vector3D InverseTransformDirection(Vector3D worldDirection)
            => Rotation.Inverse * worldDirection;

        /// <summary>
        /// Gets the inverse of this transform.
        /// </summary>
        public readonly Transform Inverse
        {
            get
            {
                var invRotation = Rotation.Inverse;
                var invScale = new Vector3D(1.0 / Scale.X, 1.0 / Scale.Y, 1.0 / Scale.Z);
                var invPosition = -(invRotation * new Vector3D(
                    Position.X * invScale.X,
                    Position.Y * invScale.Y,
                    Position.Z * invScale.Z
                ));
                return new Transform(invPosition, invRotation, invScale);
            }
        }

        /// <summary>
        /// Interpolates between two transforms.
        /// </summary>
        public static Transform Lerp(Transform a, Transform b, double t)
        {
            return new Transform(
                Vector3D.Lerp(a.Position, b.Position, t),
                Quaternion.Slerp(a.Rotation, b.Rotation, t),
                Vector3D.Lerp(a.Scale, b.Scale, t)
            );
        }

        #endregion

        #region Operators

        /// <summary>
        /// Combines two transforms.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Transform operator *(Transform parent, Transform child)
        {
            return new Transform(
                parent.TransformPoint(child.Position),
                parent.Rotation * child.Rotation,
                new Vector3D(
                    parent.Scale.X * child.Scale.X,
                    parent.Scale.Y * child.Scale.Y,
                    parent.Scale.Z * child.Scale.Z
                )
            );
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Transform a, Transform b)
            => a.Equals(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Transform a, Transform b)
            => !a.Equals(b);

        #endregion

        #region Equality

        public readonly bool Equals(Transform other)
            => Position == other.Position &&
               Rotation == other.Rotation &&
               Scale == other.Scale;

        public override readonly bool Equals(object? obj)
            => obj is Transform other && Equals(other);

        public override readonly int GetHashCode()
            => HashCode.Combine(Position, Rotation, Scale);

        public override readonly string ToString()
            => $"Transform(Pos: {Position}, Rot: {Rotation}, Scale: {Scale})";

        #endregion
    }
}
