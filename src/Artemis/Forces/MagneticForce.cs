using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Simulates magnetic attraction/repulsion between a magnetic source and magnetic bodies.
    /// Uses a simplified dipole model.
    /// </summary>
    public class MagneticForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the position of the magnet.
        /// </summary>
        public Vector3D Position { get; set; }

        /// <summary>
        /// Gets or sets the magnetic dipole moment strength (in A·m²).
        /// Positive = north pole pointing in PoleDirection.
        /// </summary>
        public double DipoleMoment { get; set; }

        /// <summary>
        /// Gets or sets the direction the magnetic pole is pointing.
        /// </summary>
        public Vector3D PoleDirection { get; set; }

        /// <summary>
        /// Gets or sets the permeability constant (μ₀/4π for scaling).
        /// Default is set for game-scale effects.
        /// </summary>
        public double Permeability { get; set; }

        /// <summary>
        /// Gets or sets the minimum distance to prevent singularity.
        /// </summary>
        public double MinDistance { get; set; }

        /// <summary>
        /// Gets or sets the maximum effective range.
        /// </summary>
        public double MaxRange { get; set; }

        /// <summary>
        /// Gets or sets whether this is an electromagnet (can be toggled).
        /// </summary>
        public bool IsElectromagnet { get; set; }

        /// <summary>
        /// Gets or sets the current for electromagnets (multiplies dipole moment).
        /// </summary>
        public double Current { get; set; }

        /// <summary>
        /// Gets or sets whether to attract (true) or repel (false).
        /// </summary>
        public bool IsAttracting { get; set; }

        /// <summary>
        /// Creates a new magnetic force.
        /// </summary>
        public MagneticForce(
            Vector3D position,
            double dipoleMoment = 100.0,
            Vector3D? poleDirection = null,
            string? id = null)
        {
            Id = id ?? $"Magnetic_{Guid.NewGuid():N}";
            Position = position;
            DipoleMoment = dipoleMoment;
            PoleDirection = poleDirection?.Normalized ?? Vector3D.Up;
            Permeability = 100.0; // Game-scale default
            MinDistance = 0.1;
            MaxRange = 20.0;
            IsElectromagnet = false;
            Current = 1.0;
            IsAttracting = true;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            var delta = position - Position;
            double distance = delta.Magnitude;

            if (distance > MaxRange || distance < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            distance = Math.Max(distance, MinDistance);

            // Calculate effective dipole moment
            double effectiveMoment = DipoleMoment;
            if (IsElectromagnet)
                effectiveMoment *= Current;

            // Simplified magnetic force: F ∝ m / r³
            // The force falls off as inverse cube of distance
            double forceMagnitude = Permeability * effectiveMoment / (distance * distance * distance);

            // Direction: towards or away from magnet
            var direction = -delta.Normalized;
            if (!IsAttracting)
                direction = -direction;

            // Add slight tangential component based on pole orientation
            var poleInfluence = Vector3D.Cross(PoleDirection, delta.Normalized);
            if (poleInfluence.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                direction = (direction + poleInfluence.Normalized * 0.1).Normalized;
            }

            return direction * forceMagnitude;
        }

        #region Fluent API

        /// <summary>
        /// Sets the magnet position.
        /// </summary>
        public MagneticForce AtPosition(Vector3D position)
        {
            Position = position;
            return this;
        }

        /// <summary>
        /// Sets the dipole moment strength.
        /// </summary>
        public MagneticForce WithStrength(double dipoleMoment)
        {
            DipoleMoment = dipoleMoment;
            return this;
        }

        /// <summary>
        /// Sets the pole direction.
        /// </summary>
        public MagneticForce WithPoleDirection(Vector3D direction)
        {
            PoleDirection = direction.Normalized;
            return this;
        }

        /// <summary>
        /// Sets the effective range.
        /// </summary>
        public MagneticForce WithRange(double maxRange)
        {
            MaxRange = maxRange;
            return this;
        }

        /// <summary>
        /// Configures as an electromagnet.
        /// </summary>
        public MagneticForce AsElectromagnet(double current = 1.0)
        {
            IsElectromagnet = true;
            Current = current;
            return this;
        }

        /// <summary>
        /// Sets attraction/repulsion mode.
        /// </summary>
        public MagneticForce Attracting(bool isAttracting = true)
        {
            IsAttracting = isAttracting;
            return this;
        }

        /// <summary>
        /// Sets repulsion mode.
        /// </summary>
        public MagneticForce Repelling()
        {
            IsAttracting = false;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a weak magnet (like a refrigerator magnet).
        /// </summary>
        public static MagneticForce Weak(Vector3D position)
            => new MagneticForce(position, 10.0).WithRange(2.0);

        /// <summary>
        /// Creates a standard magnet.
        /// </summary>
        public static MagneticForce Standard(Vector3D position)
            => new MagneticForce(position, 100.0).WithRange(10.0);

        /// <summary>
        /// Creates a powerful magnet.
        /// </summary>
        public static MagneticForce Powerful(Vector3D position)
            => new MagneticForce(position, 500.0).WithRange(20.0);

        /// <summary>
        /// Creates an industrial electromagnet.
        /// </summary>
        public static MagneticForce Industrial(Vector3D position)
            => new MagneticForce(position, 1000.0)
                .WithRange(30.0)
                .AsElectromagnet();

        #endregion
    }
}
