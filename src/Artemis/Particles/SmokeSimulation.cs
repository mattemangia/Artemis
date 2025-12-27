using System;
using Artemis.Core;
using Artemis.Forces;

namespace Artemis.Particles
{
    /// <summary>
    /// Simulates smoke using buoyant particles with turbulence and fading.
    /// </summary>
    public class SmokeSimulation
    {
        #region Fields

        private readonly ParticleSystem _particles;
        private readonly Random _random;
        private double _time;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the underlying particle system.
        /// </summary>
        public ParticleSystem Particles => _particles;

        /// <summary>
        /// Gets or sets the buoyancy strength (upward force).
        /// </summary>
        public double Buoyancy { get; set; } = 15.0;

        /// <summary>
        /// Gets or sets the turbulence strength.
        /// </summary>
        public double Turbulence { get; set; } = 5.0;

        /// <summary>
        /// Gets or sets the turbulence frequency.
        /// </summary>
        public double TurbulenceFrequency { get; set; } = 2.0;

        /// <summary>
        /// Gets or sets the horizontal spread rate.
        /// </summary>
        public double SpreadRate { get; set; } = 2.0;

        /// <summary>
        /// Gets or sets the initial particle size.
        /// </summary>
        public double InitialSize { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets the growth rate (size increase over lifetime).
        /// </summary>
        public double GrowthRate { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets the smoke color (ARGB).
        /// </summary>
        public uint SmokeColor { get; set; } = 0x80808080; // Semi-transparent gray

        /// <summary>
        /// Gets or sets the particle lifetime in seconds.
        /// </summary>
        public double Lifetime { get; set; } = 4.0;

        /// <summary>
        /// Gets or sets the emission rate (particles per second).
        /// </summary>
        public double EmissionRate { get; set; } = 50.0;

        /// <summary>
        /// Gets or sets the emission position.
        /// </summary>
        public Vector3D EmissionPosition { get; set; }

        /// <summary>
        /// Gets or sets the emission radius.
        /// </summary>
        public double EmissionRadius { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets the initial velocity.
        /// </summary>
        public Vector3D InitialVelocity { get; set; } = new(0, 2, 0);

        /// <summary>
        /// Gets or sets wind direction and strength.
        /// </summary>
        public Vector3D Wind { get; set; } = Vector3D.Zero;

        /// <summary>
        /// Gets or sets the density (affects buoyancy).
        /// Lighter smoke rises faster.
        /// </summary>
        public double Density { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets whether the emitter is active.
        /// </summary>
        public bool IsEmitting { get; set; } = true;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new smoke simulation.
        /// </summary>
        /// <param name="maxParticles">Maximum number of smoke particles.</param>
        public SmokeSimulation(int maxParticles = 5000)
        {
            _particles = new ParticleSystem(maxParticles)
            {
                Gravity = Vector3D.Zero, // Smoke ignores gravity (uses buoyancy)
                Damping = 0.98,
                DefaultLifetime = Lifetime
            };
            _random = new Random();
            _time = 0;
        }

        /// <summary>
        /// Creates a smoke simulation at a specific position.
        /// </summary>
        public SmokeSimulation(Vector3D position, int maxParticles = 5000)
            : this(maxParticles)
        {
            EmissionPosition = position;
        }

        #endregion

        #region Update

        private double _emissionAccumulator;

        /// <summary>
        /// Updates the smoke simulation.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void Update(double deltaTime)
        {
            _time += deltaTime;

            // Emit new particles
            if (IsEmitting)
            {
                _emissionAccumulator += EmissionRate * deltaTime;
                int toEmit = (int)_emissionAccumulator;
                _emissionAccumulator -= toEmit;

                for (int i = 0; i < toEmit; i++)
                {
                    EmitParticle();
                }
            }

            // Apply smoke physics to each particle
            var particles = _particles.Particles;
            for (int i = 0; i < _particles.MaxParticles; i++)
            {
                if (!particles[i].IsAlive)
                    continue;

                // Calculate age ratio (0 = just born, 1 = about to die)
                double ageRatio = 1.0 - (particles[i].Lifetime / Lifetime);
                ageRatio = Math.Clamp(ageRatio, 0, 1);

                // Buoyancy force (rises up)
                double buoyancyForce = Buoyancy * (1.0 - Density);
                particles[i].ApplyForce(new Vector3D(0, buoyancyForce * particles[i].Mass, 0));

                // Turbulence (Perlin-like noise approximation)
                double turbX = Math.Sin(_time * TurbulenceFrequency + particles[i].Position.Y * 0.5) * Turbulence;
                double turbZ = Math.Cos(_time * TurbulenceFrequency * 1.3 + particles[i].Position.Y * 0.7) * Turbulence;
                particles[i].ApplyForce(new Vector3D(turbX, 0, turbZ));

                // Horizontal spread (increases with age)
                double spreadForce = SpreadRate * ageRatio;
                Vector3D toCenter = EmissionPosition - particles[i].Position;
                toCenter.Y = 0;
                if (toCenter.MagnitudeSquared > 0.01)
                {
                    // Gentle outward push
                    particles[i].ApplyForce(-toCenter.Normalized * spreadForce);
                }

                // Wind
                if (Wind.MagnitudeSquared > 0)
                {
                    particles[i].ApplyForce(Wind * particles[i].Mass);
                }

                // Update size (smoke expands as it rises)
                particles[i].Radius = InitialSize + GrowthRate * ageRatio;

                // Update alpha (fade out)
                UpdateParticleAlpha(ref particles[i], ageRatio);
            }

            // Update particle system
            _particles.Update(deltaTime);
        }

        private void EmitParticle()
        {
            // Random position within emission radius
            double angle = _random.NextDouble() * Math.PI * 2;
            double r = _random.NextDouble() * EmissionRadius;
            var offset = new Vector3D(
                Math.Cos(angle) * r,
                0,
                Math.Sin(angle) * r
            );

            // Add some randomness to velocity
            var velocity = InitialVelocity + new Vector3D(
                (_random.NextDouble() - 0.5) * 2,
                _random.NextDouble() * 1,
                (_random.NextDouble() - 0.5) * 2
            );

            _particles.Spawn(
                EmissionPosition + offset,
                velocity,
                mass: 0.1,
                radius: InitialSize,
                lifetime: Lifetime * (0.8 + _random.NextDouble() * 0.4),
                color: SmokeColor,
                flags: ParticleFlags.None // No gravity, no world collision
            );
        }

        private void UpdateParticleAlpha(ref Particle particle, double ageRatio)
        {
            // Extract base color and calculate new alpha
            byte r = (byte)((SmokeColor >> 16) & 0xFF);
            byte g = (byte)((SmokeColor >> 8) & 0xFF);
            byte b = (byte)(SmokeColor & 0xFF);
            byte baseAlpha = (byte)((SmokeColor >> 24) & 0xFF);

            // Fade out as particle ages
            byte newAlpha = (byte)(baseAlpha * (1.0 - ageRatio * ageRatio)); // Quadratic fade

            particle.Color = (uint)((newAlpha << 24) | (r << 16) | (g << 8) | b);
        }

        #endregion

        #region Control Methods

        /// <summary>
        /// Emits a puff of smoke.
        /// </summary>
        /// <param name="count">Number of particles.</param>
        public void Puff(int count = 20)
        {
            for (int i = 0; i < count; i++)
            {
                EmitParticle();
            }
        }

        /// <summary>
        /// Clears all smoke particles.
        /// </summary>
        public void Clear()
        {
            _particles.Clear();
        }

        /// <summary>
        /// Sets the emission position.
        /// </summary>
        public SmokeSimulation At(Vector3D position)
        {
            EmissionPosition = position;
            return this;
        }

        /// <summary>
        /// Sets the wind.
        /// </summary>
        public SmokeSimulation WithWind(Vector3D wind)
        {
            Wind = wind;
            return this;
        }

        /// <summary>
        /// Sets the smoke color.
        /// </summary>
        public SmokeSimulation WithColor(byte r, byte g, byte b, byte a = 128)
        {
            SmokeColor = (uint)((a << 24) | (r << 16) | (g << 8) | b);
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a campfire smoke effect.
        /// </summary>
        public static SmokeSimulation Campfire(Vector3D position)
        {
            return new SmokeSimulation(position)
            {
                Buoyancy = 12.0,
                Turbulence = 3.0,
                EmissionRate = 30,
                InitialSize = 0.15,
                GrowthRate = 0.8,
                Lifetime = 5.0,
                SmokeColor = 0x60404040, // Dark gray, semi-transparent
                Density = 0.3
            };
        }

        /// <summary>
        /// Creates a chimney smoke effect.
        /// </summary>
        public static SmokeSimulation Chimney(Vector3D position)
        {
            return new SmokeSimulation(position)
            {
                Buoyancy = 20.0,
                Turbulence = 2.0,
                EmissionRate = 40,
                InitialSize = 0.2,
                GrowthRate = 1.2,
                Lifetime = 6.0,
                InitialVelocity = new Vector3D(0, 5, 0),
                SmokeColor = 0x50606060
            };
        }

        /// <summary>
        /// Creates a steam/vapor effect.
        /// </summary>
        public static SmokeSimulation Steam(Vector3D position)
        {
            return new SmokeSimulation(position)
            {
                Buoyancy = 8.0,
                Turbulence = 1.5,
                EmissionRate = 60,
                InitialSize = 0.08,
                GrowthRate = 0.3,
                Lifetime = 2.0,
                SmokeColor = 0x40FFFFFF, // White, very transparent
                Density = 0.1
            };
        }

        /// <summary>
        /// Creates an explosion smoke effect.
        /// </summary>
        public static SmokeSimulation Explosion(Vector3D position, int particleCount = 200)
        {
            var smoke = new SmokeSimulation(position)
            {
                Buoyancy = 25.0,
                Turbulence = 8.0,
                EmissionRate = 0, // Don't continuously emit
                InitialSize = 0.3,
                GrowthRate = 1.5,
                Lifetime = 3.0,
                SmokeColor = 0x80202020, // Dark smoke
                IsEmitting = false
            };

            // Burst of particles
            for (int i = 0; i < particleCount; i++)
            {
                var random = new Random(i);
                double theta = random.NextDouble() * Math.PI * 2;
                double phi = random.NextDouble() * Math.PI;
                double speed = 5 + random.NextDouble() * 10;

                var velocity = new Vector3D(
                    Math.Sin(phi) * Math.Cos(theta) * speed,
                    Math.Abs(Math.Cos(phi)) * speed + 3, // Bias upward
                    Math.Sin(phi) * Math.Sin(theta) * speed
                );

                smoke.Particles.Spawn(
                    position,
                    velocity,
                    mass: 0.1,
                    radius: 0.2 + random.NextDouble() * 0.3,
                    lifetime: 2.0 + random.NextDouble() * 2.0,
                    color: 0x80303030
                );
            }

            return smoke;
        }

        #endregion
    }
}
