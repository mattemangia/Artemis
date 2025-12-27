using System;
using System.Runtime.CompilerServices;
using Artemis.Core;
using Artemis.Materials;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Type of 2D physics body.
    /// </summary>
    public enum BodyType2D
    {
        /// <summary>Cannot move, infinite mass.</summary>
        Static,
        /// <summary>Moves only by explicit velocity, ignores forces.</summary>
        Kinematic,
        /// <summary>Fully simulated, responds to forces and collisions.</summary>
        Dynamic
    }

    /// <summary>
    /// Represents a rigid body in the 2D physics simulation.
    /// </summary>
    public class RigidBody2D
    {
        #region Fields

        private Vector2D _position;
        private Vector2D _velocity;
        private double _rotation;
        private double _angularVelocity;
        private double _mass;
        private double _inverseMass;
        private double _inertia;
        private double _inverseInertia;
        private Vector2D _force;
        private double _torque;
        private AABB2D _aabb;
        private double _sleepTimer;

        #endregion

        #region Properties

        /// <summary>Unique identifier for this body.</summary>
        public string Id { get; }

        /// <summary>Type of physics body.</summary>
        public BodyType2D BodyType { get; set; } = BodyType2D.Dynamic;

        /// <summary>Whether the body is active in simulation.</summary>
        public bool IsActive { get; set; } = true;

        /// <summary>Position in world space.</summary>
        public Vector2D Position
        {
            get => _position;
            set
            {
                _position = value;
                UpdateAABB();
                WakeUp();
            }
        }

        /// <summary>Linear velocity.</summary>
        public Vector2D Velocity
        {
            get => _velocity;
            set
            {
                _velocity = value;
                if (_velocity.MagnitudeSquared > PhysicsConstants.Epsilon)
                    WakeUp();
            }
        }

        /// <summary>Rotation angle in radians.</summary>
        public double Rotation
        {
            get => _rotation;
            set
            {
                _rotation = NormalizeAngle(value);
                UpdateAABB();
                WakeUp();
            }
        }

        /// <summary>Angular velocity in radians per second.</summary>
        public double AngularVelocity
        {
            get => _angularVelocity;
            set
            {
                _angularVelocity = value;
                if (Math.Abs(_angularVelocity) > PhysicsConstants.Epsilon)
                    WakeUp();
            }
        }

        /// <summary>Mass in kg.</summary>
        public double Mass
        {
            get => _mass;
            set
            {
                _mass = Math.Max(value, PhysicsConstants.Epsilon);
                _inverseMass = BodyType == BodyType2D.Dynamic ? 1.0 / _mass : 0;
            }
        }

        /// <summary>Inverse mass (0 for static/kinematic bodies).</summary>
        public double InverseMass => BodyType == BodyType2D.Dynamic ? _inverseMass : 0;

        /// <summary>Moment of inertia.</summary>
        public double Inertia
        {
            get => _inertia;
            set
            {
                _inertia = Math.Max(value, PhysicsConstants.Epsilon);
                _inverseInertia = BodyType == BodyType2D.Dynamic ? 1.0 / _inertia : 0;
            }
        }

        /// <summary>Inverse inertia (0 for static/kinematic bodies).</summary>
        public double InverseInertia => BodyType == BodyType2D.Dynamic ? _inverseInertia : 0;

        /// <summary>The collision shape attached to this body.</summary>
        public Shape2D? Shape { get; set; }

        /// <summary>Physics material (friction, restitution).</summary>
        public PhysicsMaterial Material { get; set; }

        /// <summary>Linear damping (air resistance).</summary>
        public double LinearDamping { get; set; } = 0.0;

        /// <summary>Angular damping (rotational resistance).</summary>
        public double AngularDamping { get; set; } = 0.0;

        /// <summary>Gravity scale (1 = normal, 0 = no gravity).</summary>
        public double GravityScale { get; set; } = 1.0;

        /// <summary>Whether the body is sleeping (not simulated).</summary>
        public bool IsSleeping { get; set; }

        /// <summary>Whether the body can go to sleep.</summary>
        public bool CanSleep { get; set; } = true;

        /// <summary>Whether the body is a bullet (uses continuous collision detection).</summary>
        public bool IsBullet { get; set; }

        /// <summary>Whether rotation is locked.</summary>
        public bool FixedRotation { get; set; }

        /// <summary>Axis-aligned bounding box in world space.</summary>
        public AABB2D AABB => _aabb;

        /// <summary>User-defined data.</summary>
        public object? UserData { get; set; }

        /// <summary>Collision filter category bits.</summary>
        public ushort CategoryBits { get; set; } = 0x0001;

        /// <summary>Collision filter mask bits.</summary>
        public ushort MaskBits { get; set; } = 0xFFFF;

        /// <summary>Collision group index.</summary>
        public short GroupIndex { get; set; }

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new 2D rigid body.
        /// </summary>
        public RigidBody2D(string? id = null)
        {
            Id = id ?? $"Body2D_{Guid.NewGuid():N}";
            _position = Vector2D.Zero;
            _velocity = Vector2D.Zero;
            _rotation = 0;
            _angularVelocity = 0;
            _mass = 1.0;
            _inverseMass = 1.0;
            _inertia = 1.0;
            _inverseInertia = 1.0;
            Material = new PhysicsMaterial();
            _aabb = new AABB2D(Vector2D.Zero, Vector2D.Zero);
        }

        /// <summary>
        /// Creates a new 2D rigid body at a position.
        /// </summary>
        public RigidBody2D(Vector2D position, double mass = 1.0, string? id = null)
            : this(id)
        {
            _position = position;
            Mass = mass;
        }

        #endregion

        #region Force Application

        /// <summary>
        /// Applies a force to the center of mass.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyForce(Vector2D force)
        {
            if (BodyType != BodyType2D.Dynamic)
                return;

            _force += force;
            WakeUp();
        }

        /// <summary>
        /// Applies a force at a world point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyForceAtPoint(Vector2D force, Vector2D point)
        {
            if (BodyType != BodyType2D.Dynamic)
                return;

            _force += force;
            _torque += Vector2D.Cross(point - _position, force);
            WakeUp();
        }

        /// <summary>
        /// Applies an impulse to the center of mass.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyImpulse(Vector2D impulse)
        {
            if (BodyType != BodyType2D.Dynamic)
                return;

            _velocity += impulse * _inverseMass;
            WakeUp();
        }

        /// <summary>
        /// Applies an impulse at a world point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyImpulseAtPoint(Vector2D impulse, Vector2D point)
        {
            if (BodyType != BodyType2D.Dynamic)
                return;

            _velocity += impulse * _inverseMass;
            if (!FixedRotation)
            {
                _angularVelocity += Vector2D.Cross(point - _position, impulse) * _inverseInertia;
            }
            WakeUp();
        }

        /// <summary>
        /// Applies a torque (rotational force).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyTorque(double torque)
        {
            if (BodyType != BodyType2D.Dynamic || FixedRotation)
                return;

            _torque += torque;
            WakeUp();
        }

        /// <summary>
        /// Applies an angular impulse.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyAngularImpulse(double impulse)
        {
            if (BodyType != BodyType2D.Dynamic || FixedRotation)
                return;

            _angularVelocity += impulse * _inverseInertia;
            WakeUp();
        }

        /// <summary>
        /// Clears all accumulated forces and torques.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ClearForces()
        {
            _force = Vector2D.Zero;
            _torque = 0;
        }

        #endregion

        #region Integration

        /// <summary>
        /// Integrates the body's motion for the given time step.
        /// </summary>
        public void Integrate(double dt, Vector2D gravity)
        {
            if (BodyType == BodyType2D.Static || !IsActive || IsSleeping)
                return;

            if (BodyType == BodyType2D.Dynamic)
            {
                // Apply gravity
                _velocity += gravity * GravityScale * dt;

                // Apply forces
                _velocity += _force * _inverseMass * dt;

                if (!FixedRotation)
                {
                    _angularVelocity += _torque * _inverseInertia * dt;
                }

                // Apply damping
                _velocity *= 1.0 / (1.0 + dt * LinearDamping);
                _angularVelocity *= 1.0 / (1.0 + dt * AngularDamping);
            }

            // Update position and rotation
            _position += _velocity * dt;
            _rotation = NormalizeAngle(_rotation + _angularVelocity * dt);

            // Update AABB
            UpdateAABB();

            // Check for sleep
            if (CanSleep && BodyType == BodyType2D.Dynamic)
            {
                CheckSleep(dt);
            }

            // Clear forces
            ClearForces();
        }

        /// <summary>
        /// Performs velocity integration only (for split integration).
        /// </summary>
        public void IntegrateVelocity(double dt, Vector2D gravity)
        {
            if (BodyType != BodyType2D.Dynamic || !IsActive || IsSleeping)
                return;

            // Apply gravity
            _velocity += gravity * GravityScale * dt;

            // Apply forces
            _velocity += _force * _inverseMass * dt;

            if (!FixedRotation)
            {
                _angularVelocity += _torque * _inverseInertia * dt;
            }

            // Apply damping
            _velocity *= 1.0 / (1.0 + dt * LinearDamping);
            _angularVelocity *= 1.0 / (1.0 + dt * AngularDamping);

            ClearForces();
        }

        /// <summary>
        /// Performs position integration only (for split integration).
        /// </summary>
        public void IntegratePosition(double dt)
        {
            if (BodyType == BodyType2D.Static || !IsActive || IsSleeping)
                return;

            _position += _velocity * dt;
            _rotation = NormalizeAngle(_rotation + _angularVelocity * dt);
            UpdateAABB();

            if (CanSleep && BodyType == BodyType2D.Dynamic)
            {
                CheckSleep(dt);
            }
        }

        #endregion

        #region Sleep Management

        /// <summary>
        /// Wakes up the body from sleep.
        /// </summary>
        public void WakeUp()
        {
            IsSleeping = false;
            _sleepTimer = 0;
        }

        /// <summary>
        /// Puts the body to sleep.
        /// </summary>
        public void Sleep()
        {
            if (!CanSleep)
                return;

            IsSleeping = true;
            _velocity = Vector2D.Zero;
            _angularVelocity = 0;
        }

        private void CheckSleep(double dt)
        {
            const double sleepVelocityThreshold = 0.01;
            const double sleepAngularThreshold = 0.01;
            const double sleepTimeThreshold = 0.5;

            bool isResting =
                _velocity.MagnitudeSquared < sleepVelocityThreshold * sleepVelocityThreshold &&
                Math.Abs(_angularVelocity) < sleepAngularThreshold;

            if (isResting)
            {
                _sleepTimer += dt;
                if (_sleepTimer >= sleepTimeThreshold)
                {
                    Sleep();
                }
            }
            else
            {
                _sleepTimer = 0;
            }
        }

        #endregion

        #region Utilities

        /// <summary>
        /// Updates the AABB from the shape.
        /// </summary>
        public void UpdateAABB()
        {
            if (Shape != null)
            {
                _aabb = Shape.ComputeAABB(_position, _rotation);
            }
            else
            {
                _aabb = new AABB2D(_position, _position);
            }
        }

        /// <summary>
        /// Sets the mass data from the shape and density.
        /// </summary>
        public void SetMassFromShape(double density = 1.0)
        {
            if (Shape == null)
                return;

            var massData = Shape.ComputeMass(density);
            Mass = massData.Mass;
            Inertia = massData.Inertia;
        }

        /// <summary>
        /// Transforms a world point to local coordinates.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D WorldToLocal(Vector2D worldPoint)
        {
            var d = worldPoint - _position;
            return Vector2D.Rotate(d, -_rotation);
        }

        /// <summary>
        /// Transforms a local point to world coordinates.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D LocalToWorld(Vector2D localPoint)
        {
            return _position + Vector2D.Rotate(localPoint, _rotation);
        }

        /// <summary>
        /// Gets the velocity at a world point.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D GetVelocityAtPoint(Vector2D worldPoint)
        {
            var r = worldPoint - _position;
            return _velocity + Vector2D.Cross(_angularVelocity, r);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double NormalizeAngle(double angle)
        {
            while (angle > Math.PI) angle -= 2 * Math.PI;
            while (angle < -Math.PI) angle += 2 * Math.PI;
            return angle;
        }

        /// <summary>
        /// Checks if this body should collide with another based on filters.
        /// </summary>
        public bool ShouldCollide(RigidBody2D other)
        {
            // Group filtering
            if (GroupIndex != 0 && GroupIndex == other.GroupIndex)
            {
                return GroupIndex > 0; // Positive: always collide, Negative: never collide
            }

            // Category/mask filtering
            return (CategoryBits & other.MaskBits) != 0 && (MaskBits & other.CategoryBits) != 0;
        }

        #endregion

        #region Static Factory Methods

        /// <summary>
        /// Creates a dynamic circle body.
        /// </summary>
        public static RigidBody2D CreateCircle(Vector2D position, double radius, double density = 1.0)
        {
            var body = new RigidBody2D(position)
            {
                BodyType = BodyType2D.Dynamic,
                Shape = new CircleShape(radius)
            };
            body.SetMassFromShape(density);
            body.UpdateAABB();
            return body;
        }

        /// <summary>
        /// Creates a dynamic box body.
        /// </summary>
        public static RigidBody2D CreateBox(Vector2D position, double width, double height, double density = 1.0)
        {
            var body = new RigidBody2D(position)
            {
                BodyType = BodyType2D.Dynamic,
                Shape = new BoxShape(width * 0.5, height * 0.5)
            };
            body.SetMassFromShape(density);
            body.UpdateAABB();
            return body;
        }

        /// <summary>
        /// Creates a static box body.
        /// </summary>
        public static RigidBody2D CreateStaticBox(Vector2D position, double width, double height)
        {
            var body = new RigidBody2D(position)
            {
                BodyType = BodyType2D.Static,
                Shape = new BoxShape(width * 0.5, height * 0.5)
            };
            body.UpdateAABB();
            return body;
        }

        /// <summary>
        /// Creates a static edge body.
        /// </summary>
        public static RigidBody2D CreateEdge(Vector2D start, Vector2D end)
        {
            var center = (start + end) * 0.5;
            var body = new RigidBody2D(center)
            {
                BodyType = BodyType2D.Static,
                Shape = new EdgeShape(start - center, end - center)
            };
            body.UpdateAABB();
            return body;
        }

        /// <summary>
        /// Creates a static chain body.
        /// </summary>
        public static RigidBody2D CreateChain(Vector2D[] vertices, bool loop = false)
        {
            var body = new RigidBody2D(Vector2D.Zero)
            {
                BodyType = BodyType2D.Static,
                Shape = new ChainShape(vertices, loop)
            };
            body.UpdateAABB();
            return body;
        }

        /// <summary>
        /// Creates a kinematic body (moves by velocity, unaffected by forces).
        /// </summary>
        public static RigidBody2D CreateKinematic(Vector2D position, Shape2D shape)
        {
            var body = new RigidBody2D(position)
            {
                BodyType = BodyType2D.Kinematic,
                Shape = shape
            };
            body.UpdateAABB();
            return body;
        }

        #endregion

        public override string ToString()
            => $"RigidBody2D({Id}, Pos: {Position}, Mass: {Mass:F2}kg, Type: {BodyType})";
    }
}
