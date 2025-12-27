using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Represents an aerodynamic drag force that opposes motion.
    /// Uses the drag equation: F = -0.5 * ρ * v² * Cd * A * v̂
    /// </summary>
    public class DragForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the drag coefficient (Cd).
        /// Typical values: Sphere ~0.47, Cube ~1.05, Streamlined ~0.04
        /// </summary>
        public double DragCoefficient { get; set; }

        /// <summary>
        /// Gets or sets the cross-sectional area (A) in m².
        /// </summary>
        public double CrossSectionArea { get; set; }

        /// <summary>
        /// Gets or sets the fluid density (ρ) in kg/m³.
        /// Air at sea level: ~1.225 kg/m³
        /// Water: ~1000 kg/m³
        /// </summary>
        public double FluidDensity { get; set; }

        /// <summary>
        /// Gets or sets the velocity of the fluid (for wind or current simulation).
        /// </summary>
        public Vector3D FluidVelocity { get; set; } = Vector3D.Zero;

        /// <summary>
        /// Creates a new drag force with default air properties.
        /// </summary>
        /// <param name="dragCoefficient">The drag coefficient.</param>
        /// <param name="crossSectionArea">The cross-sectional area in m².</param>
        /// <param name="id">Optional unique identifier.</param>
        public DragForce(double dragCoefficient = 0.47, double crossSectionArea = 1.0, string? id = null)
        {
            Id = id ?? $"Drag_{Guid.NewGuid():N}";
            DragCoefficient = dragCoefficient;
            CrossSectionArea = crossSectionArea;
            FluidDensity = 1.225; // Air at sea level
        }

        /// <summary>
        /// Creates a new drag force with full parameters.
        /// </summary>
        public DragForce(double dragCoefficient, double crossSectionArea, double fluidDensity, string? id = null)
        {
            Id = id ?? $"Drag_{Guid.NewGuid():N}";
            DragCoefficient = dragCoefficient;
            CrossSectionArea = crossSectionArea;
            FluidDensity = fluidDensity;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            // Relative velocity to the fluid
            var relativeVelocity = velocity - FluidVelocity;
            double speedSquared = relativeVelocity.MagnitudeSquared;

            if (speedSquared < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            double speed = Math.Sqrt(speedSquared);
            var direction = relativeVelocity / speed;

            // Drag force magnitude: F = 0.5 * ρ * v² * Cd * A
            double forceMagnitude = 0.5 * FluidDensity * speedSquared * DragCoefficient * CrossSectionArea;

            // Apply in opposite direction of motion
            return -direction * forceMagnitude;
        }

        #region Static Factory Methods

        /// <summary>
        /// Creates a drag force for air at sea level.
        /// </summary>
        public static DragForce Air(double dragCoefficient = 0.47, double crossSectionArea = 1.0, string? id = null)
            => new(dragCoefficient, crossSectionArea, 1.225, id);

        /// <summary>
        /// Creates a drag force for water.
        /// </summary>
        public static DragForce Water(double dragCoefficient = 0.47, double crossSectionArea = 1.0, string? id = null)
            => new(dragCoefficient, crossSectionArea, 1000.0, id);

        /// <summary>
        /// Creates a simple linear drag force (for simplified simulations).
        /// </summary>
        public static LinearDragForce Linear(double coefficient = 0.1, string? id = null)
            => new(coefficient, id);

        #endregion
    }

    /// <summary>
    /// Simplified linear drag force: F = -k * v
    /// Useful for game physics where realism is less important.
    /// </summary>
    public class LinearDragForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the drag coefficient.
        /// </summary>
        public double Coefficient { get; set; }

        /// <summary>
        /// Creates a new linear drag force.
        /// </summary>
        /// <param name="coefficient">The drag coefficient.</param>
        /// <param name="id">Optional unique identifier.</param>
        public LinearDragForce(double coefficient = 0.1, string? id = null)
        {
            Id = id ?? $"LinearDrag_{Guid.NewGuid():N}";
            Coefficient = coefficient;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            return -velocity * Coefficient;
        }
    }
}
