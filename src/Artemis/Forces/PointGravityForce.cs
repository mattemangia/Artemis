using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Represents gravitational attraction between masses using Newton's law of universal gravitation.
    /// F = G * (m1 * m2) / r²
    /// </summary>
    public class PointGravityForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the position of the attracting mass.
        /// </summary>
        public Vector3D Position { get; set; }

        /// <summary>
        /// Gets or sets the mass of the attracting body in kg.
        /// </summary>
        public double Mass { get; set; }

        /// <summary>
        /// Gets or sets the gravitational constant (G).
        /// Default is the real universal gravitational constant.
        /// Can be set to custom values for game physics.
        /// </summary>
        public double GravitationalConstant { get; set; }

        /// <summary>
        /// Gets or sets the minimum distance to prevent singularity at r=0.
        /// </summary>
        public double MinDistance { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets the maximum distance beyond which the force is not applied.
        /// Set to 0 or negative for infinite range.
        /// </summary>
        public double MaxDistance { get; set; } = 0;

        /// <summary>
        /// Gets or sets whether to use soft-body gravity falloff (smooth near center).
        /// </summary>
        public bool UseSoftFalloff { get; set; } = false;

        /// <summary>
        /// Creates a new point gravity force with real gravitational constant.
        /// </summary>
        /// <param name="position">Position of the attracting body.</param>
        /// <param name="mass">Mass of the attracting body in kg.</param>
        /// <param name="id">Optional unique identifier.</param>
        public PointGravityForce(Vector3D position, double mass, string? id = null)
        {
            Id = id ?? $"PointGravity_{Guid.NewGuid():N}";
            Position = position;
            Mass = mass;
            GravitationalConstant = PhysicsConstants.GravitationalConstant;
        }

        /// <summary>
        /// Creates a new point gravity force with custom gravitational constant.
        /// </summary>
        /// <param name="position">Position of the attracting body.</param>
        /// <param name="mass">Mass of the attracting body in kg.</param>
        /// <param name="gravitationalConstant">Custom gravitational constant.</param>
        /// <param name="id">Optional unique identifier.</param>
        public PointGravityForce(Vector3D position, double mass, double gravitationalConstant, string? id = null)
        {
            Id = id ?? $"PointGravity_{Guid.NewGuid():N}";
            Position = position;
            Mass = mass;
            GravitationalConstant = gravitationalConstant;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            var displacement = Position - position;
            double distanceSquared = displacement.MagnitudeSquared;
            double distance = Math.Sqrt(distanceSquared);

            // Check max distance
            if (MaxDistance > 0 && distance > MaxDistance)
                return Vector3D.Zero;

            // Prevent singularity
            if (distance < MinDistance)
            {
                if (UseSoftFalloff)
                {
                    // Smooth falloff inside minimum distance
                    double t = distance / MinDistance;
                    distanceSquared = MinDistance * MinDistance;
                    distance = MinDistance;
                    displacement = displacement.Normalized * MinDistance;

                    // Apply smooth cubic falloff
                    double smoothFactor = t * t * (3 - 2 * t);

                    var direction = displacement / distance;
                    double softForceMagnitude = GravitationalConstant * Mass * mass / distanceSquared;
                    return direction * softForceMagnitude * smoothFactor;
                }
                else
                {
                    distanceSquared = MinDistance * MinDistance;
                    distance = MinDistance;
                }
            }

            var normalizedDirection = displacement / distance;

            // Newton's law of universal gravitation: F = G * m1 * m2 / r²
            double forceMagnitude = GravitationalConstant * Mass * mass / distanceSquared;

            return normalizedDirection * forceMagnitude;
        }

        #region Static Factory Methods

        /// <summary>
        /// Creates a point gravity simulating Earth's gravity at its center.
        /// </summary>
        public static PointGravityForce Earth(Vector3D position, string? id = null)
            => new(position, 5.972e24, id);

        /// <summary>
        /// Creates a point gravity simulating the Moon.
        /// </summary>
        public static PointGravityForce Moon(Vector3D position, string? id = null)
            => new(position, 7.342e22, id);

        /// <summary>
        /// Creates a point gravity simulating the Sun.
        /// </summary>
        public static PointGravityForce Sun(Vector3D position, string? id = null)
            => new(position, 1.989e30, id);

        /// <summary>
        /// Creates a point gravity with game-friendly parameters.
        /// Uses a larger G constant for visible effects at smaller scales.
        /// </summary>
        /// <param name="position">Position of the attracting body.</param>
        /// <param name="mass">Mass of the attracting body.</param>
        /// <param name="strength">Gravity strength multiplier (1.0 = default game gravity).</param>
        /// <param name="id">Optional unique identifier.</param>
        public static PointGravityForce GameGravity(Vector3D position, double mass, double strength = 1.0, string? id = null)
        {
            // Use a much larger G for game-scale effects
            return new PointGravityForce(position, mass, 100.0 * strength, id);
        }

        #endregion
    }

    /// <summary>
    /// Represents a repulsive force (opposite of gravity).
    /// Useful for explosions, magnetic repulsion, or force fields.
    /// </summary>
    public class RepulsionForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the center of the repulsion field.
        /// </summary>
        public Vector3D Position { get; set; }

        /// <summary>
        /// Gets or sets the strength of the repulsion.
        /// </summary>
        public double Strength { get; set; }

        /// <summary>
        /// Gets or sets the falloff exponent (2 = inverse square, 1 = linear, 0 = constant).
        /// </summary>
        public double FalloffExponent { get; set; } = 2.0;

        /// <summary>
        /// Gets or sets the maximum range of the repulsion.
        /// </summary>
        public double MaxRange { get; set; } = 10.0;

        /// <summary>
        /// Gets or sets the minimum distance to prevent extreme forces.
        /// </summary>
        public double MinDistance { get; set; } = 0.1;

        /// <summary>
        /// Creates a new repulsion force.
        /// </summary>
        public RepulsionForce(Vector3D position, double strength = 100.0, double maxRange = 10.0, string? id = null)
        {
            Id = id ?? $"Repulsion_{Guid.NewGuid():N}";
            Position = position;
            Strength = strength;
            MaxRange = maxRange;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            var displacement = position - Position;
            double distance = displacement.Magnitude;

            if (distance > MaxRange || distance < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            distance = Math.Max(distance, MinDistance);

            var direction = displacement / displacement.Magnitude;

            // Calculate force with falloff
            double forceMagnitude = Strength / Math.Pow(distance, FalloffExponent);

            return direction * forceMagnitude;
        }
    }
}
