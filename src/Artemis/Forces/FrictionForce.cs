using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Simulates surface friction that slows down moving objects.
    /// Implements both static and kinetic friction.
    /// </summary>
    public class FrictionForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the surface normal (direction perpendicular to surface).
        /// </summary>
        public Vector3D SurfaceNormal { get; set; }

        /// <summary>
        /// Gets or sets the static friction coefficient.
        /// Must overcome this to start moving.
        /// </summary>
        public double StaticCoefficient { get; set; }

        /// <summary>
        /// Gets or sets the kinetic (dynamic) friction coefficient.
        /// Applied while moving.
        /// </summary>
        public double KineticCoefficient { get; set; }

        /// <summary>
        /// Gets or sets the surface height (Y coordinate for horizontal surfaces).
        /// </summary>
        public double SurfaceHeight { get; set; }

        /// <summary>
        /// Gets or sets the contact tolerance (distance to surface).
        /// </summary>
        public double ContactTolerance { get; set; }

        /// <summary>
        /// Gets or sets the velocity threshold below which static friction applies.
        /// </summary>
        public double StaticVelocityThreshold { get; set; }

        /// <summary>
        /// Creates a new friction force for a horizontal surface.
        /// </summary>
        public FrictionForce(
            double staticCoefficient = 0.5,
            double kineticCoefficient = 0.3,
            double surfaceHeight = 0,
            string? id = null)
        {
            Id = id ?? $"Friction_{Guid.NewGuid():N}";
            SurfaceNormal = Vector3D.Up;
            StaticCoefficient = staticCoefficient;
            KineticCoefficient = kineticCoefficient;
            SurfaceHeight = surfaceHeight;
            ContactTolerance = 0.1;
            StaticVelocityThreshold = 0.1;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            // Check if object is in contact with surface
            double distanceToSurface = Vector3D.Dot(position - new Vector3D(0, SurfaceHeight, 0), SurfaceNormal);
            if (distanceToSurface > ContactTolerance)
                return Vector3D.Zero;

            // Get velocity tangent to surface
            double normalVelocity = Vector3D.Dot(velocity, SurfaceNormal);
            var tangentVelocity = velocity - SurfaceNormal * normalVelocity;
            double tangentSpeed = tangentVelocity.Magnitude;

            if (tangentSpeed < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            // Calculate normal force (simplified: just use gravity)
            double normalForce = mass * PhysicsConstants.EarthGravity;

            // Determine friction type and coefficient
            double frictionCoeff = tangentSpeed < StaticVelocityThreshold
                ? StaticCoefficient
                : KineticCoefficient;

            // Friction force magnitude
            double frictionMagnitude = frictionCoeff * normalForce;

            // Friction direction: opposite to motion
            var frictionDirection = -tangentVelocity / tangentSpeed;

            // Clamp friction to not exceed the force needed to stop
            double maxFriction = mass * tangentSpeed / 0.016; // Based on 60fps
            frictionMagnitude = Math.Min(frictionMagnitude, maxFriction);

            return frictionDirection * frictionMagnitude;
        }

        #region Fluent API

        /// <summary>
        /// Sets the surface normal.
        /// </summary>
        public FrictionForce WithNormal(Vector3D normal)
        {
            SurfaceNormal = normal.Normalized;
            return this;
        }

        /// <summary>
        /// Sets both friction coefficients.
        /// </summary>
        public FrictionForce WithCoefficients(double staticCoeff, double kineticCoeff)
        {
            StaticCoefficient = staticCoeff;
            KineticCoefficient = kineticCoeff;
            return this;
        }

        /// <summary>
        /// Sets the surface height.
        /// </summary>
        public FrictionForce AtHeight(double height)
        {
            SurfaceHeight = height;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates friction for ice (very low).
        /// </summary>
        public static FrictionForce Ice(double surfaceHeight = 0)
            => new FrictionForce(0.1, 0.03, surfaceHeight);

        /// <summary>
        /// Creates friction for wood.
        /// </summary>
        public static FrictionForce Wood(double surfaceHeight = 0)
            => new FrictionForce(0.5, 0.3, surfaceHeight);

        /// <summary>
        /// Creates friction for concrete.
        /// </summary>
        public static FrictionForce Concrete(double surfaceHeight = 0)
            => new FrictionForce(0.6, 0.45, surfaceHeight);

        /// <summary>
        /// Creates friction for rubber on concrete.
        /// </summary>
        public static FrictionForce Rubber(double surfaceHeight = 0)
            => new FrictionForce(1.0, 0.8, surfaceHeight);

        /// <summary>
        /// Creates friction for metal on metal.
        /// </summary>
        public static FrictionForce Metal(double surfaceHeight = 0)
            => new FrictionForce(0.74, 0.57, surfaceHeight);

        /// <summary>
        /// Creates friction for sand/gravel.
        /// </summary>
        public static FrictionForce Sand(double surfaceHeight = 0)
            => new FrictionForce(0.6, 0.4, surfaceHeight);

        #endregion
    }
}
