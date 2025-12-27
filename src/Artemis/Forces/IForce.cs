using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Interface for force generators that can be applied to physics bodies.
    /// </summary>
    public interface IForce
    {
        /// <summary>
        /// Gets the unique identifier for this force.
        /// </summary>
        string Id { get; }

        /// <summary>
        /// Gets or sets whether this force is enabled.
        /// </summary>
        bool Enabled { get; set; }

        /// <summary>
        /// Calculates the force to apply to a body at the given position.
        /// </summary>
        /// <param name="position">The position of the body.</param>
        /// <param name="velocity">The current velocity of the body.</param>
        /// <param name="mass">The mass of the body.</param>
        /// <returns>The force vector to apply.</returns>
        Vector3D Calculate(Vector3D position, Vector3D velocity, double mass);
    }
}
