using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;
using Artemis.Forces;
using Artemis.Particles;
using Artemis.Destruction;

namespace Artemis.Modifiers
{
    /// <summary>
    /// Base interface for all physics modifiers.
    /// Modifiers are higher-level systems that combine forces, particles, and bodies.
    /// </summary>
    public interface IPhysicsModifier
    {
        /// <summary>Whether the modifier is enabled.</summary>
        bool Enabled { get; set; }

        /// <summary>Updates the modifier.</summary>
        void Update(double deltaTime);

        /// <summary>Applies the modifier to a body.</summary>
        void Apply(IPhysicsBody body, double deltaTime);

        /// <summary>Applies the modifier to a particle.</summary>
        void Apply(Particle particle, double deltaTime);
    }

    /// <summary>
    /// Wind modifier that affects particles, bodies, and can erode erodible objects.
    /// </summary>
    public class WindModifier : IPhysicsModifier
    {
        #region Fields

        private readonly Random _random = new();
        private Vector3D _currentDirection;
        private double _currentStrength;
        private double _gustTimer;
        private double _turbulencePhase;

        #endregion

        #region Properties

        /// <summary>Whether the modifier is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Base wind direction.</summary>
        public Vector3D BaseDirection { get; set; } = new(1, 0, 0);

        /// <summary>Base wind strength.</summary>
        public double BaseStrength { get; set; } = 5.0;

        /// <summary>Gust strength multiplier.</summary>
        public double GustStrength { get; set; } = 2.0;

        /// <summary>Gust frequency (gusts per second).</summary>
        public double GustFrequency { get; set; } = 0.5;

        /// <summary>Gust duration in seconds.</summary>
        public double GustDuration { get; set; } = 1.0;

        /// <summary>Turbulence intensity (0-1).</summary>
        public double Turbulence { get; set; } = 0.3;

        /// <summary>Turbulence frequency.</summary>
        public double TurbulenceFrequency { get; set; } = 2.0;

        /// <summary>Drag coefficient for bodies.</summary>
        public double BodyDragCoefficient { get; set; } = 0.47;

        /// <summary>Particle drag multiplier.</summary>
        public double ParticleDragMultiplier { get; set; } = 1.0;

        /// <summary>Current effective wind direction.</summary>
        public Vector3D CurrentDirection => _currentDirection;

        /// <summary>Current effective wind strength.</summary>
        public double CurrentStrength => _currentStrength;

        /// <summary>Optional bounds for the wind effect.</summary>
        public AABB? Bounds { get; set; }

        /// <summary>List of erodible bodies affected by this wind.</summary>
        public List<ErodibleBody> ErodibleBodies { get; } = new();

        /// <summary>Event when particles are eroded from erodible bodies.</summary>
        public event Action<List<Particle>>? OnParticlesEroded;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new wind modifier.
        /// </summary>
        public WindModifier() { }

        /// <summary>
        /// Creates a new wind modifier with direction and strength.
        /// </summary>
        public WindModifier(Vector3D direction, double strength)
        {
            BaseDirection = direction.Normalized;
            BaseStrength = strength;
        }

        #endregion

        #region Update

        /// <summary>
        /// Updates the wind simulation.
        /// </summary>
        public void Update(double deltaTime)
        {
            if (!Enabled) return;

            // Update gust
            _gustTimer += deltaTime;
            double gustPhase = _gustTimer * GustFrequency * 2 * Math.PI;
            double gustFactor = Math.Max(0, Math.Sin(gustPhase));
            gustFactor = Math.Pow(gustFactor, 3); // Sharpen gusts

            // Update turbulence
            _turbulencePhase += deltaTime * TurbulenceFrequency;
            double turbX = Math.Sin(_turbulencePhase * 1.3) * Turbulence;
            double turbY = Math.Sin(_turbulencePhase * 0.7 + 1) * Turbulence * 0.5;
            double turbZ = Math.Sin(_turbulencePhase * 1.1 + 2) * Turbulence;

            // Calculate current wind
            _currentStrength = BaseStrength * (1 + gustFactor * (GustStrength - 1));
            _currentDirection = (BaseDirection + new Vector3D(turbX, turbY, turbZ)).Normalized;

            // Erode erodible bodies
            foreach (var erodible in ErodibleBodies)
            {
                var eroded = erodible.ApplyWind(_currentDirection, _currentStrength, deltaTime);
                if (eroded.Count > 0)
                {
                    OnParticlesEroded?.Invoke(eroded);
                }
            }
        }

        /// <summary>
        /// Applies wind force to a body.
        /// </summary>
        public void Apply(IPhysicsBody body, double deltaTime)
        {
            if (!Enabled) return;
            if (Bounds.HasValue && !Bounds.Value.Contains(body.Position)) return;

            // Calculate wind force using drag equation
            // F = 0.5 * rho * v^2 * Cd * A
            var relativeVelocity = _currentDirection * _currentStrength - body.Velocity;
            double speed = relativeVelocity.Magnitude;

            if (speed > 0.01)
            {
                double area = EstimateArea(body);
                double airDensity = 1.225; // kg/m³
                double forceMagnitude = 0.5 * airDensity * speed * speed * BodyDragCoefficient * area;
                var force = relativeVelocity.Normalized * forceMagnitude;
                body.ApplyForce(force);
            }
        }

        /// <summary>
        /// Applies wind force to a particle.
        /// </summary>
        public void Apply(Particle particle, double deltaTime)
        {
            if (!Enabled) return;
            if (Bounds.HasValue && !Bounds.Value.Contains(particle.Position)) return;

            var relativeVelocity = _currentDirection * _currentStrength - particle.Velocity;
            double speed = relativeVelocity.Magnitude;

            if (speed > 0.01)
            {
                double area = particle.Size * particle.Size * Math.PI;
                double dragForce = 0.5 * 1.225 * speed * speed * 0.47 * area * ParticleDragMultiplier;
                var acceleration = relativeVelocity.Normalized * dragForce / particle.Mass;
                particle.Velocity += acceleration * deltaTime;
            }
        }

        private double EstimateArea(IPhysicsBody body)
        {
            var aabb = body.BoundingBox;
            var size = aabb.Size;
            // Estimate cross-section as average of face areas
            return (size.X * size.Y + size.Y * size.Z + size.X * size.Z) / 3;
        }

        #endregion

        #region Presets

        /// <summary>Creates a gentle breeze.</summary>
        public static WindModifier Breeze() => new()
        {
            BaseStrength = 2.0,
            GustStrength = 1.3,
            GustFrequency = 0.2,
            Turbulence = 0.1
        };

        /// <summary>Creates moderate wind.</summary>
        public static WindModifier Moderate() => new()
        {
            BaseStrength = 8.0,
            GustStrength = 1.8,
            GustFrequency = 0.4,
            Turbulence = 0.25
        };

        /// <summary>Creates strong wind.</summary>
        public static WindModifier Strong() => new()
        {
            BaseStrength = 15.0,
            GustStrength = 2.5,
            GustFrequency = 0.6,
            Turbulence = 0.4
        };

        /// <summary>Creates hurricane-force wind.</summary>
        public static WindModifier Hurricane() => new()
        {
            BaseStrength = 40.0,
            GustStrength = 3.0,
            GustFrequency = 0.8,
            Turbulence = 0.6
        };

        #endregion
    }

    /// <summary>
    /// Gravity modifier that affects all objects in a region.
    /// </summary>
    public class GravityModifier : IPhysicsModifier
    {
        /// <summary>Whether the modifier is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Gravity vector.</summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>Optional bounds for the gravity effect.</summary>
        public AABB? Bounds { get; set; }

        /// <summary>Falloff type for bounded gravity.</summary>
        public enum FalloffType { None, Linear, Quadratic }

        /// <summary>Falloff from center of bounds.</summary>
        public FalloffType Falloff { get; set; } = FalloffType.None;

        public void Update(double deltaTime) { }

        public void Apply(IPhysicsBody body, double deltaTime)
        {
            if (!Enabled) return;

            double strength = CalculateStrength(body.Position);
            if (strength > 0)
            {
                body.ApplyForce(Gravity * body.Mass * strength);
            }
        }

        public void Apply(Particle particle, double deltaTime)
        {
            if (!Enabled) return;

            double strength = CalculateStrength(particle.Position);
            if (strength > 0)
            {
                particle.Velocity += Gravity * strength * deltaTime;
            }
        }

        private double CalculateStrength(Vector3D position)
        {
            if (!Bounds.HasValue) return 1.0;
            if (!Bounds.Value.Contains(position)) return 0.0;

            if (Falloff == FalloffType.None) return 1.0;

            var center = Bounds.Value.Center;
            var halfSize = Bounds.Value.Size * 0.5;
            double maxDist = halfSize.Magnitude;
            double dist = (position - center).Magnitude;
            double t = dist / maxDist;

            return Falloff switch
            {
                FalloffType.Linear => 1 - t,
                FalloffType.Quadratic => (1 - t) * (1 - t),
                _ => 1.0
            };
        }

        /// <summary>Creates zero gravity zone.</summary>
        public static GravityModifier ZeroG(AABB bounds) => new()
        {
            Gravity = Vector3D.Zero,
            Bounds = bounds
        };

        /// <summary>Creates reversed gravity zone.</summary>
        public static GravityModifier Reversed(AABB bounds) => new()
        {
            Gravity = new Vector3D(0, 9.81, 0),
            Bounds = bounds
        };

        /// <summary>Creates low gravity zone (Moon-like).</summary>
        public static GravityModifier LowGravity(AABB bounds) => new()
        {
            Gravity = new Vector3D(0, -1.62, 0),
            Bounds = bounds
        };
    }

    /// <summary>
    /// Attractor modifier that pulls objects toward a point.
    /// </summary>
    public class AttractorModifier : IPhysicsModifier
    {
        /// <summary>Whether the modifier is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Attractor position.</summary>
        public Vector3D Position { get; set; }

        /// <summary>Attraction strength.</summary>
        public double Strength { get; set; } = 100.0;

        /// <summary>Maximum range of attraction.</summary>
        public double Range { get; set; } = 20.0;

        /// <summary>Whether to repel instead of attract.</summary>
        public bool Repel { get; set; } = false;

        /// <summary>Force falloff type.</summary>
        public enum ForceFalloff { Constant, Linear, InverseSquare }
        public ForceFalloff Falloff { get; set; } = ForceFalloff.InverseSquare;

        /// <summary>Kill zone radius (objects too close are stopped).</summary>
        public double KillRadius { get; set; } = 0.0;

        public void Update(double deltaTime) { }

        public void Apply(IPhysicsBody body, double deltaTime)
        {
            if (!Enabled) return;

            var toAttractor = Position - body.Position;
            double distance = toAttractor.Magnitude;

            if (distance < 0.001 || distance > Range) return;
            if (distance < KillRadius)
            {
                body.Velocity = Vector3D.Zero;
                return;
            }

            double forceMag = CalculateForce(distance);
            if (Repel) forceMag = -forceMag;

            var force = toAttractor.Normalized * forceMag;
            body.ApplyForce(force);
        }

        public void Apply(Particle particle, double deltaTime)
        {
            if (!Enabled) return;

            var toAttractor = Position - particle.Position;
            double distance = toAttractor.Magnitude;

            if (distance < 0.001 || distance > Range) return;

            double forceMag = CalculateForce(distance);
            if (Repel) forceMag = -forceMag;

            var acceleration = toAttractor.Normalized * forceMag / particle.Mass;
            particle.Velocity += acceleration * deltaTime;
        }

        private double CalculateForce(double distance)
        {
            return Falloff switch
            {
                ForceFalloff.Constant => Strength,
                ForceFalloff.Linear => Strength * (1 - distance / Range),
                ForceFalloff.InverseSquare => Strength / (distance * distance),
                _ => Strength
            };
        }
    }

    /// <summary>
    /// Turbulence modifier that adds random motion.
    /// </summary>
    public class TurbulenceModifier : IPhysicsModifier
    {
        private readonly Random _random = new();
        private double _time;

        /// <summary>Whether the modifier is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Turbulence strength.</summary>
        public double Strength { get; set; } = 5.0;

        /// <summary>Frequency of turbulence (oscillations per second).</summary>
        public double Frequency { get; set; } = 2.0;

        /// <summary>Spatial scale of turbulence.</summary>
        public double Scale { get; set; } = 1.0;

        /// <summary>Optional bounds.</summary>
        public AABB? Bounds { get; set; }

        public void Update(double deltaTime)
        {
            _time += deltaTime;
        }

        public void Apply(IPhysicsBody body, double deltaTime)
        {
            if (!Enabled) return;
            if (Bounds.HasValue && !Bounds.Value.Contains(body.Position)) return;

            var force = CalculateTurbulence(body.Position);
            body.ApplyForce(force * body.Mass);
        }

        public void Apply(Particle particle, double deltaTime)
        {
            if (!Enabled) return;
            if (Bounds.HasValue && !Bounds.Value.Contains(particle.Position)) return;

            var force = CalculateTurbulence(particle.Position);
            particle.Velocity += force * deltaTime;
        }

        private Vector3D CalculateTurbulence(Vector3D position)
        {
            // Simplex noise-like turbulence
            double x = position.X * Scale + _time * Frequency;
            double y = position.Y * Scale + _time * Frequency * 0.7;
            double z = position.Z * Scale + _time * Frequency * 1.3;

            return new Vector3D(
                Math.Sin(x * 1.1 + y * 0.7) * Math.Cos(z * 0.9) * Strength,
                Math.Sin(y * 0.9 + z * 1.1) * Math.Cos(x * 0.7) * Strength,
                Math.Sin(z * 0.7 + x * 0.9) * Math.Cos(y * 1.1) * Strength
            );
        }
    }

    /// <summary>
    /// Vortex/whirlwind modifier.
    /// </summary>
    public class VortexModifier : IPhysicsModifier
    {
        /// <summary>Whether the modifier is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Center position of the vortex.</summary>
        public Vector3D Position { get; set; }

        /// <summary>Axis of rotation.</summary>
        public Vector3D Axis { get; set; } = Vector3D.Up;

        /// <summary>Rotational strength.</summary>
        public double RotationalStrength { get; set; } = 20.0;

        /// <summary>Inward pull strength.</summary>
        public double InwardStrength { get; set; } = 5.0;

        /// <summary>Upward lift strength.</summary>
        public double LiftStrength { get; set; } = 10.0;

        /// <summary>Outer radius of the vortex.</summary>
        public double OuterRadius { get; set; } = 10.0;

        /// <summary>Inner radius (eye of the storm).</summary>
        public double InnerRadius { get; set; } = 1.0;

        public void Update(double deltaTime) { }

        public void Apply(IPhysicsBody body, double deltaTime)
        {
            if (!Enabled) return;

            var force = CalculateVortexForce(body.Position);
            body.ApplyForce(force);
        }

        public void Apply(Particle particle, double deltaTime)
        {
            if (!Enabled) return;

            var force = CalculateVortexForce(particle.Position);
            particle.Velocity += force / particle.Mass * deltaTime;
        }

        private Vector3D CalculateVortexForce(Vector3D position)
        {
            // Project onto plane perpendicular to axis
            var toPos = position - Position;
            var axisComponent = Vector3D.Dot(toPos, Axis) * Axis;
            var radial = toPos - axisComponent;
            double distance = radial.Magnitude;

            if (distance < 0.001 || distance > OuterRadius)
                return Vector3D.Zero;

            // Calculate strength based on distance
            double t = (distance - InnerRadius) / (OuterRadius - InnerRadius);
            t = Math.Clamp(t, 0, 1);
            double strength = Math.Sin(t * Math.PI); // Peak at middle

            // Tangential force (rotation)
            var tangent = Vector3D.Cross(Axis, radial.Normalized);
            var rotationalForce = tangent * RotationalStrength * strength;

            // Inward force
            var inwardForce = -radial.Normalized * InwardStrength * strength;

            // Lift force
            var liftForce = Axis * LiftStrength * strength;

            return rotationalForce + inwardForce + liftForce;
        }

        /// <summary>Creates a tornado.</summary>
        public static VortexModifier Tornado(Vector3D position) => new()
        {
            Position = position,
            RotationalStrength = 50.0,
            InwardStrength = 15.0,
            LiftStrength = 30.0,
            OuterRadius = 15.0,
            InnerRadius = 2.0
        };

        /// <summary>Creates a whirlpool.</summary>
        public static VortexModifier Whirlpool(Vector3D position) => new()
        {
            Position = position,
            Axis = Vector3D.Down,
            RotationalStrength = 30.0,
            InwardStrength = 20.0,
            LiftStrength = -5.0, // Pull down
            OuterRadius = 8.0,
            InnerRadius = 1.0
        };
    }

    /// <summary>
    /// Modifier that slows objects (like water or mud).
    /// </summary>
    public class DragModifier : IPhysicsModifier
    {
        /// <summary>Whether the modifier is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Drag coefficient (0-1, 1 = complete stop).</summary>
        public double DragCoefficient { get; set; } = 0.1;

        /// <summary>Optional bounds.</summary>
        public AABB? Bounds { get; set; }

        public void Update(double deltaTime) { }

        public void Apply(IPhysicsBody body, double deltaTime)
        {
            if (!Enabled) return;
            if (Bounds.HasValue && !Bounds.Value.Contains(body.Position)) return;

            // Apply drag force opposite to velocity
            var dragForce = -body.Velocity * DragCoefficient * body.Mass / deltaTime;
            body.ApplyForce(dragForce);
        }

        public void Apply(Particle particle, double deltaTime)
        {
            if (!Enabled) return;
            if (Bounds.HasValue && !Bounds.Value.Contains(particle.Position)) return;

            particle.Velocity *= (1 - DragCoefficient);
        }

        /// <summary>Creates water drag.</summary>
        public static DragModifier Water(AABB bounds) => new()
        {
            Bounds = bounds,
            DragCoefficient = 0.3
        };

        /// <summary>Creates mud/quicksand drag.</summary>
        public static DragModifier Mud(AABB bounds) => new()
        {
            Bounds = bounds,
            DragCoefficient = 0.7
        };
    }

    /// <summary>
    /// Manages multiple modifiers and applies them to physics objects.
    /// </summary>
    public class ModifierSystem
    {
        private readonly List<IPhysicsModifier> _modifiers = new();

        /// <summary>Gets all modifiers.</summary>
        public IReadOnlyList<IPhysicsModifier> Modifiers => _modifiers;

        /// <summary>Adds a modifier to the system.</summary>
        public void Add(IPhysicsModifier modifier)
        {
            _modifiers.Add(modifier);
        }

        /// <summary>Removes a modifier from the system.</summary>
        public bool Remove(IPhysicsModifier modifier)
        {
            return _modifiers.Remove(modifier);
        }

        /// <summary>Updates all modifiers.</summary>
        public void Update(double deltaTime)
        {
            foreach (var modifier in _modifiers)
            {
                modifier.Update(deltaTime);
            }
        }

        /// <summary>Applies all modifiers to a body.</summary>
        public void ApplyTo(IPhysicsBody body, double deltaTime)
        {
            foreach (var modifier in _modifiers)
            {
                modifier.Apply(body, deltaTime);
            }
        }

        /// <summary>Applies all modifiers to a particle.</summary>
        public void ApplyTo(Particle particle, double deltaTime)
        {
            foreach (var modifier in _modifiers)
            {
                modifier.Apply(particle, deltaTime);
            }
        }

        /// <summary>Applies all modifiers to multiple bodies.</summary>
        public void ApplyTo(IEnumerable<IPhysicsBody> bodies, double deltaTime)
        {
            foreach (var body in bodies)
            {
                ApplyTo(body, deltaTime);
            }
        }

        /// <summary>Applies all modifiers to a particle system.</summary>
        public void ApplyTo(ParticleSystem particleSystem, double deltaTime)
        {
            foreach (var particle in particleSystem.Particles)
            {
                ApplyTo(particle, deltaTime);
            }
        }

        /// <summary>Clears all modifiers.</summary>
        public void Clear()
        {
            _modifiers.Clear();
        }
    }
}
