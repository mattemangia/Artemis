using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Simulates an explosion force with shockwave propagation.
    /// </summary>
    public class ExplosionForce : IForce
    {
        private double _elapsedTime;
        private bool _hasExploded;

        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the explosion center.
        /// </summary>
        public Vector3D Center { get; set; }

        /// <summary>
        /// Gets or sets the initial explosion force.
        /// </summary>
        public double Force { get; set; }

        /// <summary>
        /// Gets or sets the maximum blast radius.
        /// </summary>
        public double Radius { get; set; }

        /// <summary>
        /// Gets or sets the upward bias (simulates ground reflection).
        /// 0 = spherical, 1 = fully upward.
        /// </summary>
        public double UpwardBias { get; set; }

        /// <summary>
        /// Gets or sets the duration of the explosion effect in seconds.
        /// </summary>
        public double Duration { get; set; }

        /// <summary>
        /// Gets or sets the falloff curve (1 = linear, 2 = quadratic, 3 = cubic).
        /// </summary>
        public double FalloffPower { get; set; }

        /// <summary>
        /// Gets or sets the shockwave speed (radius expansion per second).
        /// Set to 0 for instant explosion.
        /// </summary>
        public double ShockwaveSpeed { get; set; }

        /// <summary>
        /// Gets or sets the shockwave thickness.
        /// </summary>
        public double ShockwaveThickness { get; set; }

        /// <summary>
        /// Gets whether the explosion is still active.
        /// </summary>
        public bool IsActive => _elapsedTime < Duration;

        /// <summary>
        /// Gets the current shockwave radius.
        /// </summary>
        public double CurrentRadius => ShockwaveSpeed > 0
            ? Math.Min(_elapsedTime * ShockwaveSpeed, Radius)
            : Radius;

        /// <summary>
        /// Creates a new explosion force.
        /// </summary>
        public ExplosionForce(
            Vector3D center,
            double force = 1000.0,
            double radius = 10.0,
            string? id = null)
        {
            Id = id ?? $"Explosion_{Guid.NewGuid():N}";
            Center = center;
            Force = force;
            Radius = radius;
            UpwardBias = 0.3;
            Duration = 0.5;
            FalloffPower = 2.0;
            ShockwaveSpeed = 0; // Instant by default
            ShockwaveThickness = 2.0;
            _elapsedTime = 0;
            _hasExploded = false;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled || !IsActive)
                return Vector3D.Zero;

            var delta = position - Center;
            double distance = delta.Magnitude;

            if (distance > Radius || distance < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            // Handle shockwave propagation
            if (ShockwaveSpeed > 0)
            {
                double currentRadius = CurrentRadius;
                double innerRadius = Math.Max(0, currentRadius - ShockwaveThickness);

                // Only affect objects within shockwave ring
                if (distance > currentRadius || distance < innerRadius)
                    return Vector3D.Zero;
            }

            // Calculate force magnitude with falloff
            double normalizedDistance = distance / Radius;
            double falloff = Math.Pow(1.0 - normalizedDistance, FalloffPower);

            // Time decay
            double timeDecay = 1.0 - (_elapsedTime / Duration);

            double forceMagnitude = Force * falloff * timeDecay;

            // Direction: radially outward
            var direction = delta / distance;

            // Apply upward bias
            if (UpwardBias > 0)
            {
                direction = Vector3D.Lerp(direction, Vector3D.Up, UpwardBias).Normalized;
            }

            return direction * forceMagnitude;
        }

        /// <summary>
        /// Updates the explosion state. Call once per frame.
        /// </summary>
        /// <param name="deltaTime">Time since last update.</param>
        public void Update(double deltaTime)
        {
            _elapsedTime += deltaTime;
            _hasExploded = true;
        }

        /// <summary>
        /// Resets the explosion to trigger again.
        /// </summary>
        public void Reset()
        {
            _elapsedTime = 0;
            _hasExploded = false;
        }

        /// <summary>
        /// Triggers the explosion at a new position.
        /// </summary>
        public void Explode(Vector3D center)
        {
            Center = center;
            Reset();
            Enabled = true;
        }

        #region Fluent API

        /// <summary>
        /// Sets the explosion force.
        /// </summary>
        public ExplosionForce WithForce(double force)
        {
            Force = force;
            return this;
        }

        /// <summary>
        /// Sets the explosion radius.
        /// </summary>
        public ExplosionForce WithRadius(double radius)
        {
            Radius = radius;
            return this;
        }

        /// <summary>
        /// Sets the upward bias.
        /// </summary>
        public ExplosionForce WithUpwardBias(double bias)
        {
            UpwardBias = Math.Clamp(bias, 0, 1);
            return this;
        }

        /// <summary>
        /// Sets the explosion duration.
        /// </summary>
        public ExplosionForce WithDuration(double duration)
        {
            Duration = duration;
            return this;
        }

        /// <summary>
        /// Enables shockwave propagation.
        /// </summary>
        public ExplosionForce WithShockwave(double speed, double thickness = 2.0)
        {
            ShockwaveSpeed = speed;
            ShockwaveThickness = thickness;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a small explosion (grenade-like).
        /// </summary>
        public static ExplosionForce Grenade(Vector3D center)
            => new ExplosionForce(center, 500, 5)
                .WithDuration(0.3)
                .WithUpwardBias(0.4);

        /// <summary>
        /// Creates a medium explosion.
        /// </summary>
        public static ExplosionForce Medium(Vector3D center)
            => new ExplosionForce(center, 1000, 10)
                .WithDuration(0.5)
                .WithUpwardBias(0.3);

        /// <summary>
        /// Creates a large explosion with shockwave.
        /// </summary>
        public static ExplosionForce Large(Vector3D center)
            => new ExplosionForce(center, 2000, 20)
                .WithDuration(1.0)
                .WithShockwave(50, 3);

        /// <summary>
        /// Creates a nuclear-style explosion.
        /// </summary>
        public static ExplosionForce Nuclear(Vector3D center)
            => new ExplosionForce(center, 10000, 100)
                .WithDuration(3.0)
                .WithShockwave(100, 10)
                .WithUpwardBias(0.6);

        /// <summary>
        /// Creates an implosion (pulls inward).
        /// </summary>
        public static ExplosionForce Implosion(Vector3D center)
            => new ExplosionForce(center, -2000, 15)
                .WithDuration(0.8)
                .WithUpwardBias(0);

        #endregion
    }
}
