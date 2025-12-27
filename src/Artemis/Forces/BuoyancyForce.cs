using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Simulates buoyancy (Archimedes' principle) for objects in fluids.
    /// F = ρ_fluid × V_submerged × g
    /// </summary>
    public class BuoyancyForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the fluid surface height (Y coordinate).
        /// </summary>
        public double SurfaceHeight { get; set; }

        /// <summary>
        /// Gets or sets the fluid density in kg/m³.
        /// Water = 1000, Salt water = 1025, Oil = 900, Mercury = 13600
        /// </summary>
        public double FluidDensity { get; set; }

        /// <summary>
        /// Gets or sets the gravity magnitude (positive value).
        /// </summary>
        public double Gravity { get; set; }

        /// <summary>
        /// Gets or sets the object volume for buoyancy calculation.
        /// If null, volume is estimated from mass and assumed density.
        /// </summary>
        public double? ObjectVolume { get; set; }

        /// <summary>
        /// Gets or sets the assumed object density for volume estimation (kg/m³).
        /// Only used when ObjectVolume is null.
        /// </summary>
        public double AssumedObjectDensity { get; set; }

        /// <summary>
        /// Gets or sets the drag coefficient for objects moving through fluid.
        /// </summary>
        public double FluidDrag { get; set; }

        /// <summary>
        /// Gets or sets whether to apply horizontal current.
        /// </summary>
        public Vector3D CurrentVelocity { get; set; }

        /// <summary>
        /// Creates a new buoyancy force.
        /// </summary>
        public BuoyancyForce(
            double surfaceHeight = 0,
            double fluidDensity = 1000.0,
            string? id = null)
        {
            Id = id ?? $"Buoyancy_{Guid.NewGuid():N}";
            SurfaceHeight = surfaceHeight;
            FluidDensity = fluidDensity;
            Gravity = PhysicsConstants.EarthGravity;
            ObjectVolume = null;
            AssumedObjectDensity = 2000.0; // Assume denser than water by default
            FluidDrag = 0.5;
            CurrentVelocity = Vector3D.Zero;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            // Check if object is at or below surface
            if (position.Y >= SurfaceHeight)
                return Vector3D.Zero; // Above water, no buoyancy

            // Calculate submerged depth (simplified - assumes spherical object)
            double volume = ObjectVolume ?? (mass / AssumedObjectDensity);
            double radius = Math.Pow(3.0 * volume / (4.0 * Math.PI), 1.0 / 3.0);

            // Calculate submersion ratio
            double depth = SurfaceHeight - position.Y;
            double submersionRatio = Math.Min(depth / (2 * radius), 1.0);

            // Buoyancy force: F = ρ × V × g
            double submergedVolume = volume * submersionRatio;
            double buoyancyMagnitude = FluidDensity * submergedVolume * Gravity;

            var force = new Vector3D(0, buoyancyMagnitude, 0);

            // Add fluid drag when submerged
            if (FluidDrag > 0 && velocity.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                var relativeVelocity = velocity - CurrentVelocity;
                double speed = relativeVelocity.Magnitude;
                var dragDirection = -relativeVelocity.Normalized;

                // Quadratic drag: F = 0.5 × ρ × v² × Cd × A
                double crossSection = Math.PI * radius * radius * submersionRatio;
                double dragMagnitude = 0.5 * FluidDensity * speed * speed * FluidDrag * crossSection;

                force += dragDirection * dragMagnitude;
            }

            // Add current force when submerged
            if (CurrentVelocity.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                force += CurrentVelocity * (FluidDensity * 0.01 * submersionRatio);
            }

            return force;
        }

        #region Fluent API

        /// <summary>
        /// Sets the surface height.
        /// </summary>
        public BuoyancyForce AtSurface(double height)
        {
            SurfaceHeight = height;
            return this;
        }

        /// <summary>
        /// Sets the fluid density.
        /// </summary>
        public BuoyancyForce WithDensity(double density)
        {
            FluidDensity = density;
            return this;
        }

        /// <summary>
        /// Sets the fluid drag coefficient.
        /// </summary>
        public BuoyancyForce WithDrag(double drag)
        {
            FluidDrag = drag;
            return this;
        }

        /// <summary>
        /// Sets a current in the fluid.
        /// </summary>
        public BuoyancyForce WithCurrent(Vector3D current)
        {
            CurrentVelocity = current;
            return this;
        }

        /// <summary>
        /// Sets the object volume explicitly.
        /// </summary>
        public BuoyancyForce ForVolume(double volume)
        {
            ObjectVolume = volume;
            return this;
        }

        /// <summary>
        /// Sets the assumed object density for volume calculation.
        /// </summary>
        public BuoyancyForce WithObjectDensity(double density)
        {
            AssumedObjectDensity = density;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates buoyancy for fresh water.
        /// </summary>
        public static BuoyancyForce FreshWater(double surfaceHeight = 0)
            => new BuoyancyForce(surfaceHeight, 1000.0);

        /// <summary>
        /// Creates buoyancy for salt water (ocean).
        /// </summary>
        public static BuoyancyForce SaltWater(double surfaceHeight = 0)
            => new BuoyancyForce(surfaceHeight, 1025.0);

        /// <summary>
        /// Creates buoyancy for oil.
        /// </summary>
        public static BuoyancyForce Oil(double surfaceHeight = 0)
            => new BuoyancyForce(surfaceHeight, 900.0);

        /// <summary>
        /// Creates buoyancy for mercury (very dense).
        /// </summary>
        public static BuoyancyForce Mercury(double surfaceHeight = 0)
            => new BuoyancyForce(surfaceHeight, 13600.0);

        /// <summary>
        /// Creates buoyancy for air (for balloons, etc).
        /// </summary>
        public static BuoyancyForce Air(double surfaceHeight = double.MaxValue)
            => new BuoyancyForce(surfaceHeight, 1.225);

        #endregion
    }
}
