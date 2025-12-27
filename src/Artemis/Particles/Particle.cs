using System;
using Artemis.Core;
using Artemis.Materials;

namespace Artemis.Particles
{
    /// <summary>
    /// Represents a single particle in a particle system.
    /// Optimized for mass simulations with minimal overhead.
    /// </summary>
    public struct Particle
    {
        /// <summary>Position in world space.</summary>
        public Vector3D Position;

        /// <summary>Previous position (for Verlet integration).</summary>
        public Vector3D PreviousPosition;

        /// <summary>Velocity.</summary>
        public Vector3D Velocity;

        /// <summary>Accumulated force.</summary>
        public Vector3D Force;

        /// <summary>Mass of the particle.</summary>
        public double Mass;

        /// <summary>Inverse mass (1/mass). 0 for static particles.</summary>
        public double InverseMass;

        /// <summary>Radius of the particle.</summary>
        public double Radius;

        /// <summary>Lifetime remaining in seconds.</summary>
        public double Lifetime;

        /// <summary>Initial lifetime for ratio calculations.</summary>
        public double InitialLifetime;

        /// <summary>Color (packed as ARGB).</summary>
        public uint Color;

        /// <summary>Whether the particle is alive.</summary>
        public bool IsAlive;

        /// <summary>Custom user data index.</summary>
        public int UserData;

        /// <summary>Particle flags.</summary>
        public ParticleFlags Flags;

        /// <summary>
        /// Gets the lifetime ratio (0 = dead, 1 = just spawned).
        /// </summary>
        public readonly double LifetimeRatio =>
            InitialLifetime > 0 ? Lifetime / InitialLifetime : 0;

        /// <summary>
        /// Creates a new particle.
        /// </summary>
        public static Particle Create(
            Vector3D position,
            Vector3D velocity,
            double mass,
            double radius,
            double lifetime = double.PositiveInfinity,
            uint color = 0xFFFFFFFF)
        {
            return new Particle
            {
                Position = position,
                PreviousPosition = position,
                Velocity = velocity,
                Force = Vector3D.Zero,
                Mass = mass,
                InverseMass = mass > 0 ? 1.0 / mass : 0,
                Radius = radius,
                Lifetime = lifetime,
                InitialLifetime = lifetime,
                Color = color,
                IsAlive = true,
                UserData = 0,
                Flags = ParticleFlags.None
            };
        }

        /// <summary>
        /// Applies a force to the particle.
        /// </summary>
        public void ApplyForce(Vector3D force)
        {
            Force += force;
        }

        /// <summary>
        /// Applies an impulse to the particle.
        /// </summary>
        public void ApplyImpulse(Vector3D impulse)
        {
            Velocity += impulse * InverseMass;
        }

        /// <summary>
        /// Integrates the particle using semi-implicit Euler.
        /// </summary>
        public void IntegrateEuler(double deltaTime, double damping = 0.99)
        {
            if (!IsAlive || InverseMass == 0)
                return;

            // Update velocity
            Velocity += Force * InverseMass * deltaTime;
            Velocity *= damping;

            // Update position
            PreviousPosition = Position;
            Position += Velocity * deltaTime;

            // Clear forces
            Force = Vector3D.Zero;

            // Update lifetime
            if (!double.IsPositiveInfinity(Lifetime))
            {
                Lifetime -= deltaTime;
                if (Lifetime <= 0)
                    IsAlive = false;
            }
        }

        /// <summary>
        /// Integrates the particle using Verlet integration.
        /// Better for constraints and stable simulations.
        /// </summary>
        public void IntegrateVerlet(double deltaTime, double damping = 0.99)
        {
            if (!IsAlive || InverseMass == 0)
                return;

            var acceleration = Force * InverseMass;

            // Verlet integration
            var newPosition = Position * (1 + damping) - PreviousPosition * damping +
                             acceleration * deltaTime * deltaTime;

            // Update velocity (for collision response)
            Velocity = (newPosition - Position) / deltaTime;

            PreviousPosition = Position;
            Position = newPosition;

            // Clear forces
            Force = Vector3D.Zero;

            // Update lifetime
            if (!double.IsPositiveInfinity(Lifetime))
            {
                Lifetime -= deltaTime;
                if (Lifetime <= 0)
                    IsAlive = false;
            }
        }

        /// <summary>
        /// Kills the particle.
        /// </summary>
        public void Kill()
        {
            IsAlive = false;
            Lifetime = 0;
        }
    }

    /// <summary>
    /// Particle behavior flags.
    /// </summary>
    [Flags]
    public enum ParticleFlags
    {
        None = 0,
        CollidesWithWorld = 1 << 0,
        CollidesWithParticles = 1 << 1,
        AffectedByGravity = 1 << 2,
        FadesWithLifetime = 1 << 3,
        ScalesWithLifetime = 1 << 4,
        RotatesWithVelocity = 1 << 5,
        Bounces = 1 << 6,
        Sticky = 1 << 7
    }
}
