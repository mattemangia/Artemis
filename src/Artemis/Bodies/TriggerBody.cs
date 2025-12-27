using System;
using System.Collections.Generic;
using Artemis.Core;
using Artemis.Materials;

namespace Artemis.Bodies
{
    /// <summary>
    /// A body that detects overlaps but doesn't physically interact.
    /// Useful for detecting when objects enter/exit regions.
    /// </summary>
    public class TriggerBody : IPhysicsBody
    {
        private readonly HashSet<string> _overlappingBodies;
        private Vector3D _position;
        private Quaternion _rotation;
        private AABB _boundingBox;

        #region Properties

        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public BodyType BodyType
        {
            get => BodyType.Static;
            set { } // Triggers are always static
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
            }
        }

        /// <inheritdoc/>
        public Vector3D Velocity
        {
            get => Vector3D.Zero;
            set { }
        }

        /// <inheritdoc/>
        public Vector3D Acceleration
        {
            get => Vector3D.Zero;
            set { }
        }

        /// <inheritdoc/>
        public Quaternion Rotation
        {
            get => _rotation;
            set
            {
                _rotation = value;
                UpdateBoundingBox();
            }
        }

        /// <inheritdoc/>
        public Vector3D AngularVelocity
        {
            get => Vector3D.Zero;
            set { }
        }

        /// <inheritdoc/>
        public double Mass
        {
            get => 0;
            set { }
        }

        /// <inheritdoc/>
        public double InverseMass => 0;

        /// <inheritdoc/>
        public PhysicsMaterial Material { get; set; } = new();

        /// <inheritdoc/>
        public double LinearDamping { get; set; }

        /// <inheritdoc/>
        public double AngularDamping { get; set; }

        /// <inheritdoc/>
        public Transform Transform => new(_position, _rotation, Vector3D.One);

        /// <inheritdoc/>
        public AABB BoundingBox => _boundingBox;

        /// <inheritdoc/>
        public bool IsSleeping { get; set; }

        /// <inheritdoc/>
        public bool CanSleep { get; set; } = false;

        /// <inheritdoc/>
        public object? UserData { get; set; }

        /// <summary>
        /// Gets or sets the half-extents for box shape.
        /// </summary>
        public Vector3D HalfExtents { get; set; } = Vector3D.One;

        /// <summary>
        /// Gets or sets the radius for sphere shape.
        /// </summary>
        public double Radius { get; set; } = 1.0;

        /// <summary>
        /// Gets or sets the shape type.
        /// </summary>
        public CollisionShapeType ShapeType { get; set; } = CollisionShapeType.Box;

        /// <summary>
        /// Gets the IDs of currently overlapping bodies.
        /// </summary>
        public IReadOnlySet<string> OverlappingBodies => _overlappingBodies;

        /// <summary>
        /// Gets the number of overlapping bodies.
        /// </summary>
        public int OverlapCount => _overlappingBodies.Count;

        /// <summary>
        /// Gets or sets the collision filter.
        /// </summary>
        public CollisionFilter Filter { get; set; } = CollisionFilter.ForLayer(CollisionLayers.Trigger);

        #endregion

        #region Events

        /// <summary>
        /// Raised when a body enters the trigger.
        /// </summary>
        public event Action<IPhysicsBody>? OnTriggerEnter;

        /// <summary>
        /// Raised when a body exits the trigger.
        /// </summary>
        public event Action<IPhysicsBody>? OnTriggerExit;

        /// <summary>
        /// Raised every frame a body remains in the trigger.
        /// </summary>
        public event Action<IPhysicsBody>? OnTriggerStay;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new trigger body.
        /// </summary>
        public TriggerBody(string? id = null)
        {
            Id = id ?? $"Trigger_{Guid.NewGuid():N}";
            _position = Vector3D.Zero;
            _rotation = Quaternion.Identity;
            _overlappingBodies = new HashSet<string>();
            UpdateBoundingBox();
        }

        /// <summary>
        /// Creates a box trigger at the specified position.
        /// </summary>
        public TriggerBody(Vector3D position, Vector3D halfExtents, string? id = null)
            : this(id)
        {
            _position = position;
            HalfExtents = halfExtents;
            ShapeType = CollisionShapeType.Box;
            UpdateBoundingBox();
        }

        /// <summary>
        /// Creates a sphere trigger at the specified position.
        /// </summary>
        public TriggerBody(Vector3D position, double radius, string? id = null)
            : this(id)
        {
            _position = position;
            Radius = radius;
            ShapeType = CollisionShapeType.Sphere;
            UpdateBoundingBox();
        }

        #endregion

        #region Methods

        /// <summary>
        /// Checks if a body overlaps this trigger.
        /// </summary>
        public bool CheckOverlap(IPhysicsBody body)
        {
            if (!IsActive || !body.IsActive)
                return false;

            // AABB check first
            if (!_boundingBox.Intersects(body.BoundingBox))
                return false;

            // Detailed shape check
            switch (ShapeType)
            {
                case CollisionShapeType.Box:
                    return CheckBoxOverlap(body);

                case CollisionShapeType.Sphere:
                    return CheckSphereOverlap(body);

                default:
                    return _boundingBox.Intersects(body.BoundingBox);
            }
        }

        private bool CheckBoxOverlap(IPhysicsBody body)
        {
            // Simple AABB vs AABB for now
            return _boundingBox.Intersects(body.BoundingBox);
        }

        private bool CheckSphereOverlap(IPhysicsBody body)
        {
            // Get closest point on body's AABB to sphere center
            var closest = body.BoundingBox.ClosestPoint(_position);
            return Vector3D.DistanceSquared(_position, closest) <= Radius * Radius;
        }

        /// <summary>
        /// Updates overlap state for a body. Call this from physics world.
        /// </summary>
        public void UpdateOverlap(IPhysicsBody body)
        {
            bool wasOverlapping = _overlappingBodies.Contains(body.Id);
            bool isOverlapping = CheckOverlap(body);

            if (isOverlapping && !wasOverlapping)
            {
                _overlappingBodies.Add(body.Id);
                OnTriggerEnter?.Invoke(body);
            }
            else if (!isOverlapping && wasOverlapping)
            {
                _overlappingBodies.Remove(body.Id);
                OnTriggerExit?.Invoke(body);
            }
            else if (isOverlapping && wasOverlapping)
            {
                OnTriggerStay?.Invoke(body);
            }
        }

        /// <summary>
        /// Clears all overlap state.
        /// </summary>
        public void ClearOverlaps()
        {
            _overlappingBodies.Clear();
        }

        /// <inheritdoc/>
        public void ApplyForce(Vector3D force) { }

        /// <inheritdoc/>
        public void ApplyForceAtPoint(Vector3D force, Vector3D point) { }

        /// <inheritdoc/>
        public void ApplyImpulse(Vector3D impulse) { }

        /// <inheritdoc/>
        public void ApplyImpulseAtPoint(Vector3D impulse, Vector3D point) { }

        /// <inheritdoc/>
        public void ApplyTorque(Vector3D torque) { }

        /// <inheritdoc/>
        public void ClearForces() { }

        /// <inheritdoc/>
        public void Integrate(double deltaTime) { }

        /// <inheritdoc/>
        public void WakeUp() { }

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
                default:
                    _boundingBox = new AABB(_position - HalfExtents, _position + HalfExtents);
                    break;
            }
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a box trigger.
        /// </summary>
        public static TriggerBody CreateBox(Vector3D position, Vector3D halfExtents)
            => new(position, halfExtents);

        /// <summary>
        /// Creates a sphere trigger.
        /// </summary>
        public static TriggerBody CreateSphere(Vector3D position, double radius)
            => new(position, radius);

        /// <summary>
        /// Creates a trigger zone (large box for area detection).
        /// </summary>
        public static TriggerBody CreateZone(Vector3D min, Vector3D max)
        {
            var center = (min + max) * 0.5;
            var halfExtents = (max - min) * 0.5;
            return new TriggerBody(center, halfExtents);
        }

        #endregion
    }

    /// <summary>
    /// A sensor body that applies effects to overlapping bodies.
    /// Useful for force fields, damage zones, etc.
    /// </summary>
    public class SensorBody : TriggerBody
    {
        /// <summary>
        /// Gets or sets the force applied to bodies inside the sensor.
        /// </summary>
        public Vector3D Force { get; set; }

        /// <summary>
        /// Gets or sets the damping applied to bodies inside the sensor.
        /// </summary>
        public double Damping { get; set; }

        /// <summary>
        /// Gets or sets whether to scale force by overlap depth.
        /// </summary>
        public bool ScaleByDepth { get; set; }

        /// <summary>
        /// Creates a new sensor body.
        /// </summary>
        public SensorBody(Vector3D position, Vector3D halfExtents, string? id = null)
            : base(position, halfExtents, id ?? $"Sensor_{Guid.NewGuid():N}")
        {
            Filter = CollisionFilter.ForLayer(CollisionLayers.Sensor);
        }

        /// <summary>
        /// Applies sensor effects to an overlapping body.
        /// </summary>
        public void ApplyEffect(IPhysicsBody body)
        {
            if (!CheckOverlap(body))
                return;

            if (body.BodyType != BodyType.Dynamic)
                return;

            // Apply force
            if (Force.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                body.ApplyForce(Force);
            }

            // Apply damping
            if (Damping > 0)
            {
                body.Velocity *= (1 - Damping);
            }
        }

        #region Factory Methods

        /// <summary>
        /// Creates a wind zone.
        /// </summary>
        public static SensorBody WindZone(Vector3D position, Vector3D size, Vector3D windForce)
        {
            var sensor = new SensorBody(position, size * 0.5)
            {
                Force = windForce
            };
            return sensor;
        }

        /// <summary>
        /// Creates a slow zone (like water or mud).
        /// </summary>
        public static SensorBody SlowZone(Vector3D position, Vector3D size, double slowFactor = 0.1)
        {
            var sensor = new SensorBody(position, size * 0.5)
            {
                Damping = slowFactor
            };
            return sensor;
        }

        /// <summary>
        /// Creates an anti-gravity zone.
        /// </summary>
        public static SensorBody AntiGravityZone(Vector3D position, Vector3D size, double strength = 9.81)
        {
            var sensor = new SensorBody(position, size * 0.5)
            {
                Force = new Vector3D(0, strength, 0) // Upward force countering gravity
            };
            return sensor;
        }

        #endregion
    }
}
