using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Represents a uniform gravitational force field.
    /// </summary>
    public class GravityForce : IForce
    {
        private Vector3D _gravity;

        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the gravity vector (acceleration, not force).
        /// </summary>
        public Vector3D Gravity
        {
            get => _gravity;
            set => _gravity = value;
        }

        /// <summary>
        /// Gets or sets the gravity magnitude (convenience property for Y-axis gravity).
        /// Positive values point downward (negative Y).
        /// </summary>
        public double Magnitude
        {
            get => -_gravity.Y;
            set => _gravity = new Vector3D(0, -value, 0);
        }

        /// <summary>
        /// Creates a new gravity force with Earth's standard gravity.
        /// </summary>
        /// <param name="id">Optional unique identifier.</param>
        public GravityForce(string? id = null)
        {
            Id = id ?? $"Gravity_{System.Guid.NewGuid():N}";
            _gravity = new Vector3D(0, -PhysicsConstants.EarthGravity, 0);
        }

        /// <summary>
        /// Creates a new gravity force with the specified gravity vector.
        /// </summary>
        /// <param name="gravity">The gravity acceleration vector.</param>
        /// <param name="id">Optional unique identifier.</param>
        public GravityForce(Vector3D gravity, string? id = null)
        {
            Id = id ?? $"Gravity_{System.Guid.NewGuid():N}";
            _gravity = gravity;
        }

        /// <summary>
        /// Creates a new gravity force with the specified magnitude (downward).
        /// </summary>
        /// <param name="magnitude">The gravity magnitude (positive = downward).</param>
        /// <param name="id">Optional unique identifier.</param>
        public GravityForce(double magnitude, string? id = null)
        {
            Id = id ?? $"Gravity_{System.Guid.NewGuid():N}";
            _gravity = new Vector3D(0, -magnitude, 0);
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            // F = m * g
            return _gravity * mass;
        }

        #region Static Factory Methods

        /// <summary>
        /// Creates a gravity force using Earth's gravity.
        /// </summary>
        public static GravityForce Earth(string? id = null)
            => new(PhysicsConstants.EarthGravity, id);

        /// <summary>
        /// Creates a gravity force using Moon's gravity.
        /// </summary>
        public static GravityForce Moon(string? id = null)
            => new(PhysicsConstants.MoonGravity, id);

        /// <summary>
        /// Creates a gravity force using Mars' gravity.
        /// </summary>
        public static GravityForce Mars(string? id = null)
            => new(PhysicsConstants.MarsGravity, id);

        /// <summary>
        /// Creates a zero-gravity environment.
        /// </summary>
        public static GravityForce ZeroG(string? id = null)
            => new(0, id);

        #endregion
    }
}
