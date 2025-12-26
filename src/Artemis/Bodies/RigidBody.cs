using System;
using Artemis.Core;
using Artemis.Materials;

namespace Artemis.Bodies
{
    /// <summary>
    /// Represents a rigid body in the physics simulation.
    /// </summary>
    public class RigidBody : IPhysicsBody
    {
        #region Fields

        private Vector3D _position;
        private Vector3D _velocity;
        private Vector3D _acceleration;
        private Quaternion _rotation;
        private Vector3D _angularVelocity;
        private double _mass;
        private double _inverseMass;
        private Vector3D _accumulatedForce;
        private Vector3D _accumulatedTorque;
        private Matrix3x3 _inertiaTensor;
        private Matrix3x3 _inverseInertiaTensor;
        private Matrix3x3 _worldInverseInertiaTensor;
        private AABB _boundingBox;
        private double _sleepTimer;

        #endregion

        #region Properties

        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public BodyType BodyType { get; set; } = BodyType.Dynamic;

        /// <inheritdoc/>
        public bool IsActive { get; set; } = true;

        /// <inheritdoc/>
        public Vector3D Position
        {
            get => _position;
            set
            {
                _position = value;
                UpdateBoundingBox();
                WakeUp();
            }
        }

        /// <inheritdoc/>
        public Vector3D Velocity
        {
            get => _velocity;
            set
            {
                _velocity = value;
                if (_velocity.MagnitudeSquared > PhysicsConstants.Epsilon)
                    WakeUp();
            }
        }

        /// <inheritdoc/>
        public Vector3D Acceleration
        {
            get => _acceleration;
            set => _acceleration = value;
        }

        /// <inheritdoc/>
        public Quaternion Rotation
        {
            get => _rotation;
            set
            {
                _rotation = value.Normalized;
                UpdateWorldInertiaTensor();
                UpdateBoundingBox();
                WakeUp();
            }
        }

        /// <inheritdoc/>
        public Vector3D AngularVelocity
        {
            get => _angularVelocity;
            set
            {
                _angularVelocity = value;
                if (_angularVelocity.MagnitudeSquared > PhysicsConstants.Epsilon)
                    WakeUp();
            }
        }

        /// <inheritdoc/>
        public double Mass
        {
            get => _mass;
            set
            {
                _mass = Math.Max(value, PhysicsConstants.Epsilon);
                _inverseMass = BodyType == BodyType.Static ? 0 : 1.0 / _mass;
                UpdateInertiaTensor();
            }
        }

        /// <inheritdoc/>
        public double InverseMass => BodyType == BodyType.Static ? 0 : _inverseMass;

        /// <inheritdoc/>
        public PhysicsMaterial Material { get; set; }

        /// <inheritdoc/>
        public double LinearDamping { get; set; } = PhysicsConstants.DefaultLinearDamping;

        /// <inheritdoc/>
        public double AngularDamping { get; set; } = PhysicsConstants.DefaultAngularDamping;

        /// <inheritdoc/>
        public Transform Transform => new(_position, _rotation, Vector3D.One);

        /// <inheritdoc/>
        public AABB BoundingBox => _boundingBox;

        /// <inheritdoc/>
        public bool IsSleeping { get; set; }

        /// <inheritdoc/>
        public bool CanSleep { get; set; } = true;

        /// <inheritdoc/>
        public object? UserData { get; set; }

        /// <summary>
        /// Gets or sets the local inertia tensor.
        /// </summary>
        public Matrix3x3 InertiaTensor
        {
            get => _inertiaTensor;
            set
            {
                _inertiaTensor = value;
                _inverseInertiaTensor = value.Inverse;
                UpdateWorldInertiaTensor();
            }
        }

        /// <summary>
        /// Gets the inverse inertia tensor in world space.
        /// </summary>
        public Matrix3x3 WorldInverseInertiaTensor => _worldInverseInertiaTensor;

        /// <summary>
        /// Gets or sets the half-extents (for box collision shape).
        /// </summary>
        public Vector3D HalfExtents { get; set; } = new(0.5, 0.5, 0.5);

        /// <summary>
        /// Gets or sets the radius (for sphere collision shape).
        /// </summary>
        public double Radius { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets the collision shape type.
        /// </summary>
        public CollisionShapeType ShapeType { get; set; } = CollisionShapeType.Sphere;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new rigid body with default properties.
        /// </summary>
        /// <param name="id">Optional unique identifier.</param>
        public RigidBody(string? id = null)
        {
            Id = id ?? $"RigidBody_{Guid.NewGuid():N}";
            _position = Vector3D.Zero;
            _velocity = Vector3D.Zero;
            _acceleration = Vector3D.Zero;
            _rotation = Quaternion.Identity;
            _angularVelocity = Vector3D.Zero;
            _mass = 1.0;
            _inverseMass = 1.0;
            Material = new PhysicsMaterial();
            _inertiaTensor = Matrix3x3.Identity;
            _inverseInertiaTensor = Matrix3x3.Identity;
            _worldInverseInertiaTensor = Matrix3x3.Identity;
            UpdateBoundingBox();
        }

        /// <summary>
        /// Creates a new rigid body at the specified position.
        /// </summary>
        /// <param name="position">Initial position.</param>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="id">Optional unique identifier.</param>
        public RigidBody(Vector3D position, double mass = 1.0, string? id = null)
            : this(id)
        {
            _position = position;
            Mass = mass;
            UpdateBoundingBox();
        }

        #endregion

        #region Methods

        /// <inheritdoc/>
        public void ApplyForce(Vector3D force)
        {
            if (BodyType == BodyType.Static)
                return;

            _accumulatedForce += force;
            WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyForceAtPoint(Vector3D force, Vector3D point)
        {
            if (BodyType == BodyType.Static)
                return;

            _accumulatedForce += force;

            // Calculate torque: τ = r × F
            var r = point - _position;
            _accumulatedTorque += Vector3D.Cross(r, force);
            WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyImpulse(Vector3D impulse)
        {
            if (BodyType == BodyType.Static)
                return;

            _velocity += impulse * _inverseMass;
            WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyImpulseAtPoint(Vector3D impulse, Vector3D point)
        {
            if (BodyType == BodyType.Static)
                return;

            _velocity += impulse * _inverseMass;

            var r = point - _position;
            var angularImpulse = Vector3D.Cross(r, impulse);
            _angularVelocity += _worldInverseInertiaTensor * angularImpulse;
            WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyTorque(Vector3D torque)
        {
            if (BodyType == BodyType.Static)
                return;

            _accumulatedTorque += torque;
            WakeUp();
        }

        /// <inheritdoc/>
        public void ClearForces()
        {
            _accumulatedForce = Vector3D.Zero;
            _accumulatedTorque = Vector3D.Zero;
        }

        /// <inheritdoc/>
        public void Integrate(double deltaTime)
        {
            if (BodyType == BodyType.Static || !IsActive || IsSleeping)
                return;

            // Calculate acceleration from accumulated forces
            _acceleration = _accumulatedForce * _inverseMass;

            // Update linear velocity (semi-implicit Euler)
            _velocity += _acceleration * deltaTime;

            // Apply linear damping
            _velocity *= Math.Pow(1.0 - LinearDamping, deltaTime);

            // Update position
            _position += _velocity * deltaTime;

            // Calculate angular acceleration
            var angularAcceleration = _worldInverseInertiaTensor * _accumulatedTorque;

            // Update angular velocity
            _angularVelocity += angularAcceleration * deltaTime;

            // Apply angular damping
            _angularVelocity *= Math.Pow(1.0 - AngularDamping, deltaTime);

            // Update rotation
            if (_angularVelocity.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                var angularDelta = _angularVelocity * deltaTime;
                double angle = angularDelta.Magnitude;
                if (angle > PhysicsConstants.Epsilon)
                {
                    var axis = angularDelta / angle;
                    var deltaRotation = Quaternion.FromAxisAngle(axis, angle);
                    _rotation = (deltaRotation * _rotation).Normalized;
                }
            }

            // Update world-space inertia tensor
            UpdateWorldInertiaTensor();

            // Update bounding box
            UpdateBoundingBox();

            // Check for sleep
            CheckSleep(deltaTime);

            // Clear forces for next frame
            ClearForces();
        }

        /// <inheritdoc/>
        public void WakeUp()
        {
            IsSleeping = false;
            _sleepTimer = 0;
        }

        /// <inheritdoc/>
        public void UpdateBoundingBox()
        {
            switch (ShapeType)
            {
                case CollisionShapeType.Sphere:
                    var r = new Vector3D(Radius);
                    _boundingBox = new AABB(_position - r, _position + r);
                    break;

                case CollisionShapeType.Box:
                    // Rotate corners and find AABB
                    var corners = new Vector3D[8];
                    for (int i = 0; i < 8; i++)
                    {
                        var sign = new Vector3D(
                            (i & 1) == 0 ? -1 : 1,
                            (i & 2) == 0 ? -1 : 1,
                            (i & 4) == 0 ? -1 : 1
                        );
                        var localCorner = new Vector3D(
                            sign.X * HalfExtents.X,
                            sign.Y * HalfExtents.Y,
                            sign.Z * HalfExtents.Z
                        );
                        corners[i] = _position + _rotation * localCorner;
                    }
                    _boundingBox = AABB.FromPoints(corners);
                    break;

                default:
                    _boundingBox = new AABB(_position - Vector3D.One, _position + Vector3D.One);
                    break;
            }
        }

        private void UpdateInertiaTensor()
        {
            switch (ShapeType)
            {
                case CollisionShapeType.Sphere:
                    double sphereI = 0.4 * _mass * Radius * Radius;
                    _inertiaTensor = Matrix3x3.Diagonal(sphereI);
                    break;

                case CollisionShapeType.Box:
                    double boxIx = _mass * (HalfExtents.Y * HalfExtents.Y + HalfExtents.Z * HalfExtents.Z) / 3.0;
                    double boxIy = _mass * (HalfExtents.X * HalfExtents.X + HalfExtents.Z * HalfExtents.Z) / 3.0;
                    double boxIz = _mass * (HalfExtents.X * HalfExtents.X + HalfExtents.Y * HalfExtents.Y) / 3.0;
                    _inertiaTensor = Matrix3x3.Diagonal(boxIx, boxIy, boxIz);
                    break;

                default:
                    _inertiaTensor = Matrix3x3.Diagonal(_mass);
                    break;
            }

            _inverseInertiaTensor = _inertiaTensor.Inverse;
            UpdateWorldInertiaTensor();
        }

        private void UpdateWorldInertiaTensor()
        {
            if (BodyType == BodyType.Static)
            {
                _worldInverseInertiaTensor = Matrix3x3.Zero;
                return;
            }

            var rotMatrix = Matrix3x3.FromQuaternion(_rotation);
            _worldInverseInertiaTensor = rotMatrix * _inverseInertiaTensor * rotMatrix.Transposed;
        }

        private void CheckSleep(double deltaTime)
        {
            if (!CanSleep)
                return;

            bool isResting =
                _velocity.MagnitudeSquared < PhysicsConstants.DefaultSleepVelocityThreshold *
                                             PhysicsConstants.DefaultSleepVelocityThreshold &&
                _angularVelocity.MagnitudeSquared < PhysicsConstants.DefaultSleepAngularVelocityThreshold *
                                                    PhysicsConstants.DefaultSleepAngularVelocityThreshold;

            if (isResting)
            {
                _sleepTimer += deltaTime;
                if (_sleepTimer >= PhysicsConstants.DefaultSleepTimeThreshold)
                {
                    IsSleeping = true;
                    _velocity = Vector3D.Zero;
                    _angularVelocity = Vector3D.Zero;
                }
            }
            else
            {
                _sleepTimer = 0;
            }
        }

        #endregion

        #region Static Factory Methods

        /// <summary>
        /// Creates a sphere rigid body.
        /// </summary>
        public static RigidBody CreateSphere(Vector3D position, double radius, double mass, PhysicsMaterial? material = null)
        {
            var body = new RigidBody(position, mass)
            {
                ShapeType = CollisionShapeType.Sphere,
                Radius = radius,
                Material = material ?? new PhysicsMaterial()
            };
            body.UpdateBoundingBox();
            return body;
        }

        /// <summary>
        /// Creates a box rigid body.
        /// </summary>
        public static RigidBody CreateBox(Vector3D position, Vector3D halfExtents, double mass, PhysicsMaterial? material = null)
        {
            var body = new RigidBody(position, mass)
            {
                ShapeType = CollisionShapeType.Box,
                HalfExtents = halfExtents,
                Material = material ?? new PhysicsMaterial()
            };
            body.UpdateBoundingBox();
            return body;
        }

        /// <summary>
        /// Creates a static box (e.g., for floors, walls).
        /// </summary>
        public static RigidBody CreateStaticBox(Vector3D position, Vector3D halfExtents, PhysicsMaterial? material = null)
        {
            var body = CreateBox(position, halfExtents, 0, material);
            body.BodyType = BodyType.Static;
            return body;
        }

        /// <summary>
        /// Creates a static sphere.
        /// </summary>
        public static RigidBody CreateStaticSphere(Vector3D position, double radius, PhysicsMaterial? material = null)
        {
            var body = CreateSphere(position, radius, 0, material);
            body.BodyType = BodyType.Static;
            return body;
        }

        #endregion

        public override string ToString()
            => $"RigidBody({Id}, Pos: {Position}, Mass: {Mass:F2}kg)";
    }

    /// <summary>
    /// Defines the collision shape type.
    /// </summary>
    public enum CollisionShapeType
    {
        Sphere,
        Box,
        Capsule,
        Cylinder,
        ConvexHull,
        TriangleMesh
    }
}
