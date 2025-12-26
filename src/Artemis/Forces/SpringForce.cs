using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Represents a spring force using Hooke's law: F = -k * (x - x₀)
    /// Can be used for elastic connections, soft body physics, or cloth simulation.
    /// </summary>
    public class SpringForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the anchor point of the spring.
        /// </summary>
        public Vector3D Anchor { get; set; }

        /// <summary>
        /// Gets or sets the spring stiffness coefficient (k) in N/m.
        /// Higher values = stiffer spring.
        /// </summary>
        public double Stiffness { get; set; }

        /// <summary>
        /// Gets or sets the rest length of the spring in meters.
        /// </summary>
        public double RestLength { get; set; }

        /// <summary>
        /// Gets or sets the damping coefficient.
        /// Reduces oscillation over time.
        /// </summary>
        public double Damping { get; set; }

        /// <summary>
        /// Gets or sets the maximum force the spring can exert.
        /// Set to 0 or negative for unlimited force.
        /// </summary>
        public double MaxForce { get; set; } = 0;

        /// <summary>
        /// Creates a new spring force.
        /// </summary>
        /// <param name="anchor">The anchor point of the spring.</param>
        /// <param name="stiffness">Spring stiffness (k) in N/m.</param>
        /// <param name="restLength">Rest length of the spring.</param>
        /// <param name="damping">Damping coefficient (optional).</param>
        /// <param name="id">Optional unique identifier.</param>
        public SpringForce(
            Vector3D anchor,
            double stiffness = 100.0,
            double restLength = 1.0,
            double damping = 0.1,
            string? id = null)
        {
            Id = id ?? $"Spring_{Guid.NewGuid():N}";
            Anchor = anchor;
            Stiffness = stiffness;
            RestLength = restLength;
            Damping = damping;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            var displacement = position - Anchor;
            double distance = displacement.Magnitude;

            if (distance < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            var direction = displacement / distance;

            // Hooke's law: F = -k * (x - x₀)
            double stretch = distance - RestLength;
            double springForceMagnitude = -Stiffness * stretch;

            // Damping force: F = -c * v (projected onto spring direction)
            double dampingForceMagnitude = -Damping * Vector3D.Dot(velocity, direction);

            double totalMagnitude = springForceMagnitude + dampingForceMagnitude;

            // Apply max force limit if set
            if (MaxForce > 0)
                totalMagnitude = Math.Clamp(totalMagnitude, -MaxForce, MaxForce);

            return direction * totalMagnitude;
        }
    }

    /// <summary>
    /// Represents a spring connection between two physics bodies.
    /// </summary>
    public class BungeeForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the anchor point.
        /// </summary>
        public Vector3D Anchor { get; set; }

        /// <summary>
        /// Gets or sets the spring stiffness.
        /// </summary>
        public double Stiffness { get; set; }

        /// <summary>
        /// Gets or sets the rest length.
        /// </summary>
        public double RestLength { get; set; }

        /// <summary>
        /// Creates a new bungee force (only pulls when stretched beyond rest length).
        /// </summary>
        public BungeeForce(Vector3D anchor, double stiffness = 100.0, double restLength = 1.0, string? id = null)
        {
            Id = id ?? $"Bungee_{Guid.NewGuid():N}";
            Anchor = anchor;
            Stiffness = stiffness;
            RestLength = restLength;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            var displacement = position - Anchor;
            double distance = displacement.Magnitude;

            // Only apply force when stretched beyond rest length
            if (distance <= RestLength)
                return Vector3D.Zero;

            var direction = displacement / distance;
            double stretch = distance - RestLength;
            double forceMagnitude = -Stiffness * stretch;

            return direction * forceMagnitude;
        }
    }
}
