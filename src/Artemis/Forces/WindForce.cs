using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Simulates wind force with turbulence and gusts.
    /// </summary>
    public class WindForce : IForce
    {
        private readonly Random _random;
        private double _gustTimer;
        private double _currentGustStrength;
        private Vector3D _currentGustDirection;

        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the base wind direction (normalized automatically).
        /// </summary>
        public Vector3D Direction { get; set; }

        /// <summary>
        /// Gets or sets the base wind strength in N.
        /// </summary>
        public double Strength { get; set; }

        /// <summary>
        /// Gets or sets the turbulence factor (0 = smooth, 1 = very turbulent).
        /// </summary>
        public double Turbulence { get; set; }

        /// <summary>
        /// Gets or sets the gust frequency (gusts per second).
        /// </summary>
        public double GustFrequency { get; set; }

        /// <summary>
        /// Gets or sets the maximum gust strength multiplier.
        /// </summary>
        public double GustStrength { get; set; }

        /// <summary>
        /// Gets or sets the drag coefficient affecting how much velocity reduces wind effect.
        /// </summary>
        public double DragCoefficient { get; set; }

        /// <summary>
        /// Gets or sets whether the force scales with cross-sectional area (mass as proxy).
        /// </summary>
        public bool ScaleWithMass { get; set; }

        /// <summary>
        /// Creates a new wind force.
        /// </summary>
        public WindForce(
            Vector3D direction,
            double strength = 10.0,
            double turbulence = 0.2,
            string? id = null)
        {
            Id = id ?? $"Wind_{Guid.NewGuid():N}";
            Direction = direction.Normalized;
            Strength = strength;
            Turbulence = turbulence;
            GustFrequency = 0.5;
            GustStrength = 2.0;
            DragCoefficient = 0.1;
            ScaleWithMass = false;
            _random = new Random();
            _gustTimer = 0;
            _currentGustStrength = 0;
            _currentGustDirection = direction.Normalized;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            // Base wind direction
            var windDir = Direction.Normalized;

            // Add turbulence
            if (Turbulence > 0)
            {
                var turbulenceOffset = new Vector3D(
                    (_random.NextDouble() - 0.5) * Turbulence,
                    (_random.NextDouble() - 0.5) * Turbulence * 0.5,
                    (_random.NextDouble() - 0.5) * Turbulence
                );
                windDir = (windDir + turbulenceOffset).Normalized;
            }

            // Calculate effective wind strength
            double effectiveStrength = Strength * (1 + _currentGustStrength);

            // Reduce effect based on velocity in wind direction
            double relativeSpeed = Vector3D.Dot(velocity, windDir);
            if (relativeSpeed > 0 && DragCoefficient > 0)
            {
                effectiveStrength *= Math.Max(0, 1 - relativeSpeed * DragCoefficient);
            }

            // Scale with mass if enabled
            if (ScaleWithMass)
            {
                effectiveStrength *= Math.Pow(mass, 0.333); // Cube root for volume-based scaling
            }

            return windDir * effectiveStrength;
        }

        /// <summary>
        /// Updates gust simulation (call once per frame).
        /// </summary>
        /// <param name="deltaTime">Time since last update.</param>
        public void UpdateGusts(double deltaTime)
        {
            _gustTimer += deltaTime;

            if (_gustTimer >= 1.0 / GustFrequency && GustFrequency > 0)
            {
                _gustTimer = 0;
                _currentGustStrength = _random.NextDouble() * GustStrength;

                // Random gust direction variation
                var gustOffset = new Vector3D(
                    (_random.NextDouble() - 0.5) * 0.3,
                    (_random.NextDouble() - 0.5) * 0.1,
                    (_random.NextDouble() - 0.5) * 0.3
                );
                _currentGustDirection = (Direction + gustOffset).Normalized;
            }
            else
            {
                // Decay gust
                _currentGustStrength *= Math.Pow(0.1, deltaTime);
            }
        }

        #region Fluent API

        /// <summary>
        /// Sets the wind direction.
        /// </summary>
        public WindForce WithDirection(Vector3D direction)
        {
            Direction = direction.Normalized;
            return this;
        }

        /// <summary>
        /// Sets the wind strength.
        /// </summary>
        public WindForce WithStrength(double strength)
        {
            Strength = strength;
            return this;
        }

        /// <summary>
        /// Sets the turbulence level.
        /// </summary>
        public WindForce WithTurbulence(double turbulence)
        {
            Turbulence = Math.Clamp(turbulence, 0, 1);
            return this;
        }

        /// <summary>
        /// Configures gust behavior.
        /// </summary>
        public WindForce WithGusts(double frequency, double strength)
        {
            GustFrequency = frequency;
            GustStrength = strength;
            return this;
        }

        /// <summary>
        /// Sets whether force scales with mass.
        /// </summary>
        public WindForce WithMassScaling(bool enabled = true)
        {
            ScaleWithMass = enabled;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a gentle breeze.
        /// </summary>
        public static WindForce Breeze(Vector3D direction)
            => new WindForce(direction, 5.0, 0.1)
                .WithGusts(0.2, 0.5);

        /// <summary>
        /// Creates moderate wind.
        /// </summary>
        public static WindForce Moderate(Vector3D direction)
            => new WindForce(direction, 20.0, 0.3)
                .WithGusts(0.5, 1.5);

        /// <summary>
        /// Creates strong wind.
        /// </summary>
        public static WindForce Strong(Vector3D direction)
            => new WindForce(direction, 50.0, 0.5)
                .WithGusts(1.0, 2.5);

        /// <summary>
        /// Creates storm-level wind.
        /// </summary>
        public static WindForce Storm(Vector3D direction)
            => new WindForce(direction, 100.0, 0.8)
                .WithGusts(2.0, 4.0);

        #endregion
    }
}
