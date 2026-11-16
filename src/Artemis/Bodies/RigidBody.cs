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
        private BodyType _bodyType = BodyType.Dynamic;
        private CollisionShapeType _shapeType = CollisionShapeType.Sphere;
        private Vector3D _halfExtents = new(0.5, 0.5, 0.5);
        private double _radius = 0.5;
        private bool _hasCustomInertia;

        #endregion

        #region Properties

        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public BodyType BodyType
        {
            get => _bodyType;
            set
            {
                if (_bodyType == value)
                    return;

                _bodyType = value;
                if (value == BodyType.Static)
                {
                    _velocity = Vector3D.Zero;
                    _angularVelocity = Vector3D.Zero;
                }

                // Only dynamic bodies respond to impulses, so the world-space
                // inverse inertia has to follow the body type.
                UpdateWorldInertiaTensor();
                WakeUp();
            }
        }

        /// <summary>
        /// Gets or sets whether the body is static (alias for BodyType).
        /// </summary>
        public bool IsStatic
        {
            get => BodyType == BodyType.Static;
            set => BodyType = value ? BodyType.Static : BodyType.Dynamic;
        }

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
                if (IsSleeping)
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
                if (IsSleeping && _velocity.MagnitudeSquared > PhysicsConstants.Epsilon)
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
                if (IsSleeping)
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
                if (IsSleeping && _angularVelocity.MagnitudeSquared > PhysicsConstants.Epsilon)
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
                _inverseMass = 1.0 / _mass;
                UpdateInertiaTensor();
            }
        }

        /// <inheritdoc/>
        public double InverseMass => _bodyType == BodyType.Dynamic ? _inverseMass : 0;

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
                _hasCustomInertia = true;
                UpdateWorldInertiaTensor();
            }
        }

        /// <summary>
        /// Gets the inverse inertia tensor in world space.
        /// </summary>
        public Matrix3x3 WorldInverseInertiaTensor => _worldInverseInertiaTensor;

        /// <summary>
        /// Gets the inverse inertia tensor in world space (alias).
        /// </summary>
        public Matrix3x3 InverseInertiaTensorWorld => _worldInverseInertiaTensor;

        /// <summary>
        /// Gets the angular acceleration derived from accumulated torque.
        /// </summary>
        public Vector3D AngularAcceleration => _worldInverseInertiaTensor * _accumulatedTorque;

        /// <summary>
        /// Gets or sets the sleep timer (seconds).
        /// </summary>
        public double SleepTimer
        {
            get => _sleepTimer;
            set => _sleepTimer = value;
        }

        /// <summary>
        /// Gets or sets the half-extents (for box collision shape).
        /// </summary>
        public Vector3D HalfExtents
        {
            get => _halfExtents;
            set
            {
                _halfExtents = value;
                OnShapeChanged();
            }
        }

        /// <summary>
        /// Gets or sets the radius (for sphere collision shape).
        /// </summary>
        public double Radius
        {
            get => _radius;
            set
            {
                _radius = value;
                OnShapeChanged();
            }
        }

        /// <summary>
        /// Gets or sets the collision shape type.
        /// </summary>
        public CollisionShapeType ShapeType
        {
            get => _shapeType;
            set
            {
                _shapeType = value;
                OnShapeChanged();
            }
        }

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
            UpdateInertiaTensor();
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
            if (_bodyType != BodyType.Dynamic)
                return;

            _accumulatedForce += force;
            if (IsSleeping)
                WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyForceAtPoint(Vector3D force, Vector3D point)
        {
            if (_bodyType != BodyType.Dynamic)
                return;

            _accumulatedForce += force;

            // Calculate torque: τ = r × F
            var r = point - _position;
            _accumulatedTorque += Vector3D.Cross(r, force);
            if (IsSleeping)
                WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyImpulse(Vector3D impulse)
        {
            if (_bodyType != BodyType.Dynamic)
                return;

            _velocity += impulse * _inverseMass;
            WakeUp();
        }

        /// <inheritdoc/>
        public void ApplyImpulseAtPoint(Vector3D impulse, Vector3D point)
        {
            if (_bodyType != BodyType.Dynamic)
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
            if (_bodyType != BodyType.Dynamic)
                return;

            _accumulatedTorque += torque;
            if (IsSleeping)
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
            if (_bodyType == BodyType.Static || !IsActive || IsSleeping)
                return;

            IntegrateVelocity(deltaTime);
            IntegratePosition(deltaTime);
            CheckSleep(deltaTime);
        }

        /// <summary>
        /// First half of a simulation step: turns the accumulated forces and torques
        /// into velocity changes (semi-implicit Euler), applies damping and clears the
        /// accumulators. Positions are left untouched so that a contact solver can
        /// correct the velocities before they are integrated.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void IntegrateVelocity(double deltaTime)
        {
            if (!IsActive || IsSleeping)
                return;

            if (_bodyType != BodyType.Dynamic)
            {
                // Static and kinematic bodies are not driven by forces.
                _acceleration = Vector3D.Zero;
                ClearForces();
                return;
            }

            _acceleration = _accumulatedForce * _inverseMass;
            _velocity += _acceleration * deltaTime;
            _velocity *= Math.Pow(1.0 - LinearDamping, deltaTime);

            _angularVelocity += (_worldInverseInertiaTensor * _accumulatedTorque) * deltaTime;
            _angularVelocity *= Math.Pow(1.0 - AngularDamping, deltaTime);

            ClearForces();
        }

        /// <summary>
        /// Second half of a simulation step: advances position and orientation using
        /// the current velocities.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void IntegratePosition(double deltaTime)
        {
            if (_bodyType == BodyType.Static || !IsActive || IsSleeping)
                return;

            _position += _velocity * deltaTime;

            if (_angularVelocity.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                if (_bodyType == BodyType.Dynamic)
                    IntegrateFreeRotation(deltaTime);
                else
                    Rotate(_angularVelocity * deltaTime);

                UpdateWorldInertiaTensor();
            }

            UpdateBoundingBox();
        }

        /// <summary>
        /// Torque-free rotation of a rigid body. The angular velocity is advanced in
        /// the body frame with Euler's equations, I dw/dt = (I w) x w, using RK4. This
        /// keeps the gyroscopic coupling between the axes, so asymmetric bodies tumble
        /// and precess correctly while angular momentum and energy stay (very nearly)
        /// constant - a plain "w stays constant in world space" update conserves neither.
        /// </summary>
        private void IntegrateFreeRotation(double deltaTime)
        {
            var toLocal = _rotation.Conjugate;
            var w = toLocal * _angularVelocity;

            // Keep the rotation per sub-step small enough for RK4 to stay accurate.
            double anglePerStep = w.Magnitude * deltaTime;
            int subSteps = Math.Min(16, Math.Max(1, (int)Math.Ceiling(anglePerStep / 0.25)));
            double h = deltaTime / subSteps;

            for (int i = 0; i < subSteps; i++)
            {
                var k1 = EulerRate(w);
                var k2 = EulerRate(w + k1 * (0.5 * h));
                var k3 = EulerRate(w + k2 * (0.5 * h));
                var k4 = EulerRate(w + k3 * h);
                var next = w + (k1 + (k2 + k3) * 2.0 + k4) * (h / 6.0);

                // Orientation follows the mean body-frame angular velocity of the sub-step.
                var localDelta = (w + next) * (0.5 * h);
                double angle = localDelta.Magnitude;
                if (angle > PhysicsConstants.Epsilon)
                    _rotation = (_rotation * Quaternion.FromAxisAngle(localDelta / angle, angle)).Normalized;

                w = next;
            }

            _angularVelocity = _rotation * w;
        }

        private Vector3D EulerRate(Vector3D w)
            => _inverseInertiaTensor * Vector3D.Cross(_inertiaTensor * w, w);

        private void Rotate(Vector3D angularDelta)
        {
            double angle = angularDelta.Magnitude;
            if (angle > PhysicsConstants.Epsilon)
                _rotation = (Quaternion.FromAxisAngle(angularDelta / angle, angle) * _rotation).Normalized;
        }

        /// <summary>
        /// Writes solver results back without touching the sleep state.
        /// </summary>
        internal void SetVelocities(Vector3D linear, Vector3D angular)
        {
            _velocity = linear;
            _angularVelocity = angular;
        }

        /// <summary>
        /// Moves and turns the body without waking it or resetting its sleep timer.
        /// </summary>
        /// <param name="translation">World-space offset.</param>
        /// <param name="rotation">World-space rotation vector (axis times angle, radians).</param>
        internal void ApplyPositionCorrection(Vector3D translation, Vector3D rotation)
        {
            _position += translation;

            if (rotation.MagnitudeSquared > PhysicsConstants.Epsilon * PhysicsConstants.Epsilon)
            {
                Rotate(rotation);
                UpdateWorldInertiaTensor();
            }

            UpdateBoundingBox();
        }

        /// <inheritdoc/>
        public void WakeUp()
        {
            IsSleeping = false;
            _sleepTimer = 0;
        }

        /// <summary>
        /// Puts the body to sleep and clears velocities.
        /// </summary>
        public void Sleep()
        {
            IsSleeping = true;
            _velocity = Vector3D.Zero;
            _angularVelocity = Vector3D.Zero;
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

        private void OnShapeChanged()
        {
            UpdateInertiaTensor();
            UpdateBoundingBox();
        }

        private void UpdateInertiaTensor()
        {
            // An explicitly assigned tensor is kept until the caller replaces it.
            if (_hasCustomInertia)
                return;

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
            if (_bodyType != BodyType.Dynamic)
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
