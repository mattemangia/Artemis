using Artemis.Core;
using Artemis.Materials;

namespace Artemis.Bodies
{
    /// <summary>
    /// Defines the body type for physics simulation.
    /// </summary>
    public enum BodyType
    {
        /// <summary>
        /// A dynamic body affected by forces and collisions.
        /// </summary>
        Dynamic,

        /// <summary>
        /// A static body that doesn't move but affects other bodies.
        /// </summary>
        Static,

        /// <summary>
        /// A kinematic body controlled by user code, affects other bodies but isn't affected by forces.
        /// </summary>
        Kinematic
    }

    /// <summary>
    /// Interface for all physics bodies in the simulation.
    /// </summary>
    public interface IPhysicsBody
    {
        /// <summary>
        /// Gets the unique identifier for this body.
        /// </summary>
        string Id { get; }

        /// <summary>
        /// Gets or sets the body type.
        /// </summary>
        BodyType BodyType { get; set; }

        /// <summary>
        /// Gets or sets whether this body is active in the simulation.
        /// </summary>
        bool IsActive { get; set; }

        /// <summary>
        /// Gets or sets the position of the body's center of mass.
        /// </summary>
        Vector3D Position { get; set; }

        /// <summary>
        /// Gets or sets the linear velocity.
        /// </summary>
        Vector3D Velocity { get; set; }

        /// <summary>
        /// Gets or sets the linear acceleration.
        /// </summary>
        Vector3D Acceleration { get; set; }

        /// <summary>
        /// Gets or sets the rotation as a quaternion.
        /// </summary>
        Quaternion Rotation { get; set; }

        /// <summary>
        /// Gets or sets the angular velocity (radians per second).
        /// </summary>
        Vector3D AngularVelocity { get; set; }

        /// <summary>
        /// Gets or sets the mass in kg.
        /// </summary>
        double Mass { get; set; }

        /// <summary>
        /// Gets the inverse mass (1/mass). Zero for static bodies.
        /// </summary>
        double InverseMass { get; }

        /// <summary>
        /// Gets or sets the physical material.
        /// </summary>
        PhysicsMaterial Material { get; set; }

        /// <summary>
        /// Gets or sets the linear damping coefficient.
        /// </summary>
        double LinearDamping { get; set; }

        /// <summary>
        /// Gets or sets the angular damping coefficient.
        /// </summary>
        double AngularDamping { get; set; }

        /// <summary>
        /// Gets the transform of this body.
        /// </summary>
        Transform Transform { get; }

        /// <summary>
        /// Gets the axis-aligned bounding box of this body.
        /// </summary>
        AABB BoundingBox { get; }

        /// <summary>
        /// Gets or sets whether this body is sleeping (optimization for resting bodies).
        /// </summary>
        bool IsSleeping { get; set; }

        /// <summary>
        /// Gets or sets whether this body can sleep.
        /// </summary>
        bool CanSleep { get; set; }

        /// <summary>
        /// Gets or sets custom user data attached to this body.
        /// </summary>
        object? UserData { get; set; }

        /// <summary>
        /// Applies a force at the center of mass.
        /// </summary>
        /// <param name="force">The force to apply in Newtons.</param>
        void ApplyForce(Vector3D force);

        /// <summary>
        /// Applies a force at a specific world point.
        /// </summary>
        /// <param name="force">The force to apply in Newtons.</param>
        /// <param name="point">The point of application in world coordinates.</param>
        void ApplyForceAtPoint(Vector3D force, Vector3D point);

        /// <summary>
        /// Applies an impulse at the center of mass.
        /// </summary>
        /// <param name="impulse">The impulse to apply in N·s.</param>
        void ApplyImpulse(Vector3D impulse);

        /// <summary>
        /// Applies an impulse at a specific world point.
        /// </summary>
        /// <param name="impulse">The impulse to apply in N·s.</param>
        /// <param name="point">The point of application in world coordinates.</param>
        void ApplyImpulseAtPoint(Vector3D impulse, Vector3D point);

        /// <summary>
        /// Applies a torque.
        /// </summary>
        /// <param name="torque">The torque to apply in N·m.</param>
        void ApplyTorque(Vector3D torque);

        /// <summary>
        /// Clears all accumulated forces and torques.
        /// </summary>
        void ClearForces();

        /// <summary>
        /// Updates the body's state for the given time step.
        /// </summary>
        /// <param name="deltaTime">The time step in seconds.</param>
        void Integrate(double deltaTime);

        /// <summary>
        /// Wakes up the body from sleep.
        /// </summary>
        void WakeUp();

        /// <summary>
        /// Updates the bounding box after a transform change.
        /// </summary>
        void UpdateBoundingBox();
    }
}
