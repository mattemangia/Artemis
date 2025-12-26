using System;
using Artemis.Core;

namespace Artemis.Particles
{
    /// <summary>
    /// Simulates fire using temperature-based particles with combustion effects.
    /// </summary>
    public class FireSimulation
    {
        #region Fields

        private readonly ParticleSystem _particles;
        private readonly Random _random;
        private double _time;
        private double _emissionAccumulator;

        // Temperature tracking per particle (parallel array)
        private double[] _temperatures;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the underlying particle system.
        /// </summary>
        public ParticleSystem Particles => _particles;

        /// <summary>
        /// Gets or sets the emission position.
        /// </summary>
        public Vector3D EmissionPosition { get; set; }

        /// <summary>
        /// Gets or sets the emission radius.
        /// </summary>
        public double EmissionRadius { get; set; } = 0.3;

        /// <summary>
        /// Gets or sets the emission rate (particles per second).
        /// </summary>
        public double EmissionRate { get; set; } = 100.0;

        /// <summary>
        /// Gets or sets the particle lifetime in seconds.
        /// </summary>
        public double Lifetime { get; set; } = 1.5;

        /// <summary>
        /// Gets or sets the initial temperature in Kelvin.
        /// </summary>
        public double InitialTemperature { get; set; } = 1500.0;

        /// <summary>
        /// Gets or sets the cooling rate (temperature loss per second).
        /// </summary>
        public double CoolingRate { get; set; } = 800.0;

        /// <summary>
        /// Gets or sets the buoyancy strength.
        /// </summary>
        public double Buoyancy { get; set; } = 20.0;

        /// <summary>
        /// Gets or sets the turbulence strength.
        /// </summary>
        public double Turbulence { get; set; } = 8.0;

        /// <summary>
        /// Gets or sets the flicker frequency.
        /// </summary>
        public double FlickerFrequency { get; set; } = 10.0;

        /// <summary>
        /// Gets or sets the initial velocity.
        /// </summary>
        public Vector3D InitialVelocity { get; set; } = new(0, 3, 0);

        /// <summary>
        /// Gets or sets the velocity randomness.
        /// </summary>
        public double VelocityRandomness { get; set; } = 1.5;

        /// <summary>
        /// Gets or sets the fire intensity (0-1).
        /// </summary>
        public double Intensity { get; set; } = 1.0;

        /// <summary>
        /// Gets or sets the flame height multiplier.
        /// </summary>
        public double FlameHeight { get; set; } = 1.0;

        /// <summary>
        /// Gets or sets whether the fire is burning.
        /// </summary>
        public bool IsBurning { get; set; } = true;

        /// <summary>
        /// Gets or sets the wind affecting the fire.
        /// </summary>
        public Vector3D Wind { get; set; } = Vector3D.Zero;

        /// <summary>
        /// Gets or sets the initial particle size.
        /// </summary>
        public double InitialSize { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets whether to emit smoke at the top.
        /// </summary>
        public bool EmitSmoke { get; set; } = true;

        /// <summary>
        /// Gets the linked smoke simulation (if any).
        /// </summary>
        public SmokeSimulation? Smoke { get; private set; }

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new fire simulation.
        /// </summary>
        /// <param name="maxParticles">Maximum number of fire particles.</param>
        public FireSimulation(int maxParticles = 3000)
        {
            _particles = new ParticleSystem(maxParticles)
            {
                Gravity = Vector3D.Zero, // Fire uses buoyancy
                Damping = 0.95,
                DefaultLifetime = Lifetime
            };
            _temperatures = new double[maxParticles];
            _random = new Random();
            _time = 0;
        }

        /// <summary>
        /// Creates a fire simulation at a specific position.
        /// </summary>
        public FireSimulation(Vector3D position, int maxParticles = 3000)
            : this(maxParticles)
        {
            EmissionPosition = position;
        }

        #endregion

        #region Update

        /// <summary>
        /// Updates the fire simulation.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void Update(double deltaTime)
        {
            _time += deltaTime;

            // Emit new particles
            if (IsBurning && Intensity > 0)
            {
                _emissionAccumulator += EmissionRate * Intensity * deltaTime;
                int toEmit = (int)_emissionAccumulator;
                _emissionAccumulator -= toEmit;

                for (int i = 0; i < toEmit; i++)
                {
                    EmitParticle();
                }
            }

            // Apply fire physics to each particle
            var particles = _particles.Particles;
            for (int i = 0; i < _particles.MaxParticles; i++)
            {
                if (!particles[i].IsAlive)
                    continue;

                // Cool down the particle
                _temperatures[i] -= CoolingRate * deltaTime;

                // Kill if too cold
                if (_temperatures[i] < 400) // Below visible red
                {
                    // Optionally spawn smoke particle here
                    if (EmitSmoke && Smoke != null && _random.NextDouble() < 0.3)
                    {
                        Smoke.Particles.Spawn(
                            particles[i].Position,
                            particles[i].Velocity * 0.5,
                            mass: 0.1,
                            radius: particles[i].Radius * 1.5,
                            lifetime: 3.0,
                            color: 0x40404040
                        );
                    }
                    particles[i].Kill();
                    continue;
                }

                // Update color based on temperature
                particles[i].Color = TemperatureToColor(_temperatures[i]);

                // Buoyancy (hot air rises)
                double tempRatio = _temperatures[i] / InitialTemperature;
                double buoyancyForce = Buoyancy * tempRatio * FlameHeight;
                particles[i].ApplyForce(new Vector3D(0, buoyancyForce * particles[i].Mass, 0));

                // Turbulence with flicker
                double flicker = 1.0 + 0.3 * Math.Sin(_time * FlickerFrequency + i * 0.1);
                double turbX = Math.Sin(_time * 5 + particles[i].Position.Y * 2) * Turbulence * flicker;
                double turbZ = Math.Cos(_time * 4.3 + particles[i].Position.Y * 1.7) * Turbulence * flicker;
                particles[i].ApplyForce(new Vector3D(turbX, 0, turbZ));

                // Wind
                if (Wind.MagnitudeSquared > 0)
                {
                    particles[i].ApplyForce(Wind * particles[i].Mass);
                }

                // Shrink as it cools
                particles[i].Radius = InitialSize * tempRatio;
            }

            // Update particle system
            _particles.Update(deltaTime);

            // Update smoke if linked
            Smoke?.Update(deltaTime);
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

            // Random velocity with upward bias
            var velocity = InitialVelocity * FlameHeight + new Vector3D(
                (_random.NextDouble() - 0.5) * VelocityRandomness * 2,
                _random.NextDouble() * VelocityRandomness,
                (_random.NextDouble() - 0.5) * VelocityRandomness * 2
            );

            int index = _particles.Spawn(
                EmissionPosition + offset,
                velocity,
                mass: 0.05,
                radius: InitialSize,
                lifetime: Lifetime * (0.8 + _random.NextDouble() * 0.4),
                color: TemperatureToColor(InitialTemperature),
                flags: ParticleFlags.None
            );

            if (index >= 0)
            {
                _temperatures[index] = InitialTemperature * (0.9 + _random.NextDouble() * 0.2);
            }
        }

        #endregion

        #region Color Calculation

        /// <summary>
        /// Converts temperature to fire color using blackbody radiation approximation.
        /// </summary>
        private static uint TemperatureToColor(double temperature)
        {
            // Blackbody radiation color approximation
            // Temperature in Kelvin

            double t = temperature / 100.0;
            double r, g, b;

            // Red
            if (t <= 66)
                r = 255;
            else
            {
                r = t - 60;
                r = 329.698727446 * Math.Pow(r, -0.1332047592);
                r = Math.Clamp(r, 0, 255);
            }

            // Green
            if (t <= 66)
            {
                g = t;
                g = 99.4708025861 * Math.Log(g) - 161.1195681661;
                g = Math.Clamp(g, 0, 255);
            }
            else
            {
                g = t - 60;
                g = 288.1221695283 * Math.Pow(g, -0.0755148492);
                g = Math.Clamp(g, 0, 255);
            }

            // Blue
            if (t >= 66)
                b = 255;
            else if (t <= 19)
                b = 0;
            else
            {
                b = t - 10;
                b = 138.5177312231 * Math.Log(b) - 305.0447927307;
                b = Math.Clamp(b, 0, 255);
            }

            // Alpha based on temperature (hotter = more opaque)
            double alpha = Math.Clamp((temperature - 400) / 1000.0, 0.3, 1.0);

            return (uint)(
                ((byte)(alpha * 255) << 24) |
                ((byte)r << 16) |
                ((byte)g << 8) |
                (byte)b
            );
        }

        /// <summary>
        /// Gets a color at a specific temperature.
        /// </summary>
        public static (byte r, byte g, byte b, byte a) GetFireColor(double temperature)
        {
            uint color = TemperatureToColor(temperature);
            return (
                (byte)((color >> 16) & 0xFF),
                (byte)((color >> 8) & 0xFF),
                (byte)(color & 0xFF),
                (byte)((color >> 24) & 0xFF)
            );
        }

        #endregion

        #region Control Methods

        /// <summary>
        /// Links a smoke simulation to this fire.
        /// </summary>
        public FireSimulation WithSmoke(SmokeSimulation smoke)
        {
            Smoke = smoke;
            Smoke.EmissionPosition = EmissionPosition + new Vector3D(0, FlameHeight * 2, 0);
            Smoke.IsEmitting = false; // Fire controls smoke emission
            return this;
        }

        /// <summary>
        /// Creates and links a default smoke simulation.
        /// </summary>
        public FireSimulation WithSmoke()
        {
            var smoke = new SmokeSimulation(EmissionPosition + new Vector3D(0, FlameHeight * 2, 0))
            {
                IsEmitting = false,
                Buoyancy = 10,
                EmissionRate = 0
            };
            return WithSmoke(smoke);
        }

        /// <summary>
        /// Clears all fire particles.
        /// </summary>
        public void Extinguish()
        {
            IsBurning = false;
            _particles.Clear();
        }

        /// <summary>
        /// Ignites the fire.
        /// </summary>
        public void Ignite()
        {
            IsBurning = true;
        }

        /// <summary>
        /// Sets the fire intensity.
        /// </summary>
        public FireSimulation WithIntensity(double intensity)
        {
            Intensity = Math.Clamp(intensity, 0, 2);
            return this;
        }

        /// <summary>
        /// Sets the wind.
        /// </summary>
        public FireSimulation WithWind(Vector3D wind)
        {
            Wind = wind;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a campfire.
        /// </summary>
        public static FireSimulation Campfire(Vector3D position)
        {
            return new FireSimulation(position)
            {
                EmissionRadius = 0.4,
                EmissionRate = 80,
                InitialTemperature = 1200,
                CoolingRate = 600,
                Buoyancy = 15,
                Turbulence = 6,
                FlameHeight = 1.2,
                InitialSize = 0.12,
                Lifetime = 1.2
            }.WithSmoke();
        }

        /// <summary>
        /// Creates a torch fire.
        /// </summary>
        public static FireSimulation Torch(Vector3D position)
        {
            return new FireSimulation(position)
            {
                EmissionRadius = 0.08,
                EmissionRate = 50,
                InitialTemperature = 1400,
                CoolingRate = 900,
                Buoyancy = 25,
                Turbulence = 4,
                FlameHeight = 0.8,
                InitialSize = 0.06,
                Lifetime = 0.8
            };
        }

        /// <summary>
        /// Creates a candle flame.
        /// </summary>
        public static FireSimulation Candle(Vector3D position)
        {
            return new FireSimulation(position)
            {
                EmissionRadius = 0.02,
                EmissionRate = 30,
                InitialTemperature = 1100,
                CoolingRate = 500,
                Buoyancy = 10,
                Turbulence = 1.5,
                FlameHeight = 0.4,
                InitialSize = 0.03,
                Lifetime = 0.5,
                VelocityRandomness = 0.3
            };
        }

        /// <summary>
        /// Creates a bonfire (large fire).
        /// </summary>
        public static FireSimulation Bonfire(Vector3D position)
        {
            return new FireSimulation(position, 5000)
            {
                EmissionRadius = 1.0,
                EmissionRate = 200,
                InitialTemperature = 1600,
                CoolingRate = 500,
                Buoyancy = 30,
                Turbulence = 10,
                FlameHeight = 2.5,
                InitialSize = 0.2,
                Lifetime = 2.0
            }.WithSmoke();
        }

        /// <summary>
        /// Creates an inferno (intense fire).
        /// </summary>
        public static FireSimulation Inferno(Vector3D position)
        {
            return new FireSimulation(position, 8000)
            {
                EmissionRadius = 2.0,
                EmissionRate = 400,
                InitialTemperature = 2000,
                CoolingRate = 600,
                Buoyancy = 50,
                Turbulence = 15,
                FlameHeight = 4.0,
                InitialSize = 0.25,
                Lifetime = 2.5,
                Intensity = 1.5
            }.WithSmoke();
        }

        /// <summary>
        /// Creates a gas burner flame (blue fire).
        /// </summary>
        public static FireSimulation GasBurner(Vector3D position)
        {
            return new FireSimulation(position)
            {
                EmissionRadius = 0.05,
                EmissionRate = 60,
                InitialTemperature = 2500, // Hot enough to be blue
                CoolingRate = 1500,
                Buoyancy = 15,
                Turbulence = 2,
                FlameHeight = 0.6,
                InitialSize = 0.04,
                Lifetime = 0.4,
                EmitSmoke = false
            };
        }

        /// <summary>
        /// Creates an explosion fireball.
        /// </summary>
        public static FireSimulation Explosion(Vector3D position, double radius = 3.0)
        {
            var fire = new FireSimulation(position, 3000)
            {
                EmissionRadius = radius,
                EmissionRate = 0, // Burst only
                InitialTemperature = 2000,
                CoolingRate = 1000,
                Buoyancy = 40,
                Turbulence = 20,
                FlameHeight = 1.0,
                InitialSize = 0.3,
                Lifetime = 1.5,
                IsBurning = false
            };

            // Burst of fire particles
            var random = new Random();
            for (int i = 0; i < 500; i++)
            {
                double theta = random.NextDouble() * Math.PI * 2;
                double phi = random.NextDouble() * Math.PI;
                double speed = 5 + random.NextDouble() * 15;
                double r = random.NextDouble() * radius * 0.5;

                var offset = new Vector3D(
                    Math.Sin(phi) * Math.Cos(theta) * r,
                    Math.Abs(random.NextDouble()) * r,
                    Math.Sin(phi) * Math.Sin(theta) * r
                );

                var velocity = new Vector3D(
                    Math.Sin(phi) * Math.Cos(theta) * speed,
                    Math.Abs(Math.Cos(phi)) * speed + 5,
                    Math.Sin(phi) * Math.Sin(theta) * speed
                );

                int index = fire.Particles.Spawn(
                    position + offset,
                    velocity,
                    mass: 0.1,
                    radius: 0.2 + random.NextDouble() * 0.3,
                    lifetime: 0.8 + random.NextDouble() * 1.0,
                    color: TemperatureToColor(1500 + random.NextDouble() * 1000)
                );

                if (index >= 0)
                {
                    fire._temperatures[index] = 1500 + random.NextDouble() * 1000;
                }
            }

            return fire;
        }

        #endregion
    }
}
