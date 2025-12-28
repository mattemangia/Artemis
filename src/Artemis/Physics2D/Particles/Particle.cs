using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Artemis.Physics2D.Particles
{
    /// <summary>
    /// Types of particles in the collider.
    /// Based on real particle physics.
    /// </summary>
    public enum ParticleType
    {
        Proton,      // Massive, positively charged
        Electron,    // Light, negatively charged
        Positron,    // Antimatter electron
        Neutron,     // Neutral, heavy
        Photon,      // Massless, light speed
        Higgs,       // Heavy boson (discovered at LHC!)
        Quark,       // Fundamental particle
        Muon,        // Heavy electron
        Neutrino,    // Nearly massless, weakly interacting
        Pion         // Composite particle
    }

    /// <summary>
    /// Particle properties based on real physics.
    /// Masses in atomic mass units (amu), scaled for visualization.
    /// </summary>
    public class ParticleProperties
    {
        public ParticleType Type { get; set; }
        public double Mass { get; set; }           // Relative mass
        public int Charge { get; set; }            // Electric charge (-1, 0, +1)
        public uint Color { get; set; }            // RGBA color
        public string Symbol { get; set; } = "";
        public double Energy { get; set; }         // Kinetic energy (GeV scaled)
        public bool IsStable { get; set; }         // Does it decay?
        public double HalfLife { get; set; }       // Decay time (seconds)

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ParticleProperties Create(ParticleType type)
        {
            return type switch
            {
                ParticleType.Proton => new ParticleProperties
                {
                    Type = ParticleType.Proton,
                    Mass = 1.0,
                    Charge = +1,
                    Color = 0xFF0000FF, // Red
                    Symbol = "p⁺",
                    IsStable = true,
                    HalfLife = double.PositiveInfinity
                },
                ParticleType.Electron => new ParticleProperties
                {
                    Type = ParticleType.Electron,
                    Mass = 0.0005, // ~1/2000 of proton
                    Charge = -1,
                    Color = 0x0000FFFF, // Blue
                    Symbol = "e⁻",
                    IsStable = true,
                    HalfLife = double.PositiveInfinity
                },
                ParticleType.Positron => new ParticleProperties
                {
                    Type = ParticleType.Positron,
                    Mass = 0.0005,
                    Charge = +1,
                    Color = 0x00FFFFFF, // Cyan
                    Symbol = "e⁺",
                    IsStable = true,
                    HalfLife = double.PositiveInfinity
                },
                ParticleType.Neutron => new ParticleProperties
                {
                    Type = ParticleType.Neutron,
                    Mass = 1.0,
                    Charge = 0,
                    Color = 0x808080FF, // Gray
                    Symbol = "n⁰",
                    IsStable = false,
                    HalfLife = 881.5 // seconds
                },
                ParticleType.Photon => new ParticleProperties
                {
                    Type = ParticleType.Photon,
                    Mass = 0.0001, // Nearly massless
                    Charge = 0,
                    Color = 0xFFFF00FF, // Yellow
                    Symbol = "γ",
                    IsStable = true,
                    HalfLife = double.PositiveInfinity
                },
                ParticleType.Higgs => new ParticleProperties
                {
                    Type = ParticleType.Higgs,
                    Mass = 125.0, // GeV/c² (scaled)
                    Charge = 0,
                    Color = 0xFF00FFFF, // Magenta
                    Symbol = "H⁰",
                    IsStable = false,
                    HalfLife = 1.6e-22 // Extremely short!
                },
                ParticleType.Quark => new ParticleProperties
                {
                    Type = ParticleType.Quark,
                    Mass = 0.3,
                    Charge = -1,
                    Color = 0x00FF00FF, // Green
                    Symbol = "q",
                    IsStable = false,
                    HalfLife = 1e-23
                },
                ParticleType.Muon => new ParticleProperties
                {
                    Type = ParticleType.Muon,
                    Mass = 0.11,
                    Charge = -1,
                    Color = 0xFF8000FF, // Orange
                    Symbol = "μ⁻",
                    IsStable = false,
                    HalfLife = 2.2e-6
                },
                ParticleType.Neutrino => new ParticleProperties
                {
                    Type = ParticleType.Neutrino,
                    Mass = 0.00001,
                    Charge = 0,
                    Color = 0xFFFFFFFF, // White
                    Symbol = "ν",
                    IsStable = true,
                    HalfLife = double.PositiveInfinity
                },
                ParticleType.Pion => new ParticleProperties
                {
                    Type = ParticleType.Pion,
                    Mass = 0.14,
                    Charge = 0,
                    Color = 0x8000FFFF, // Purple
                    Symbol = "π⁰",
                    IsStable = false,
                    HalfLife = 8.4e-17
                },
                _ => throw new ArgumentException($"Unknown particle type: {type}")
            };
        }
    }

    /// <summary>
    /// Wrapper around RigidBody2D with particle-specific properties.
    /// Implements ICharged for magnetic field interactions.
    /// </summary>
    public class Particle2D : ICharged
    {
        private static readonly Random Rng = new();
        private static int _nextId;

        public int Id { get; }
        public RigidBody2D Body { get; }
        public ParticleProperties Properties { get; }
        public double Age { get; private set; }
        public bool HasDecayed { get; private set; }
        public List<Particle2D> DecayProducts { get; } = new();
        public Vector2D InitialPosition { get; }
        public double TotalEnergyDeposited { get; set; }

        // Trail for visualization
        public List<Vector2D> Trail { get; } = new();
        public int MaxTrailLength { get; set; } = 50;

        // ICharged implementation
        public double Charge => Properties.Charge;

        public Particle2D(ParticleType type, Vector2D position, Vector2D velocity)
        {
            Id = Interlocked.Increment(ref _nextId);
            Properties = ParticleProperties.Create(type);

            // Create physics body with small radius for realistic particle appearance
            double visualRadius = Math.Max(0.15, Properties.Mass * 0.2);
            Body = RigidBody2D.CreateCircle(position, visualRadius, Properties.Mass / (Math.PI * visualRadius * visualRadius));
            Body.Velocity = velocity;
            Body.Material = new Materials.PhysicsMaterial { Restitution = 1.0, Friction = 0.0 }; // Elastic, no friction
            Body.IsBullet = true;   // Use CCD for high-speed particles
            Body.UserData = this;   // Link back to particle
            Body.LinearDamping = 0; // No damping in vacuum
            Body.CanSleep = false;  // Particles should never sleep

            InitialPosition = position;

            // Calculate initial kinetic energy
            Properties.Energy = 0.5 * Properties.Mass * velocity.MagnitudeSquared;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Update(double deltaTime)
        {
            Age += deltaTime;

            // Record trail
            if (Trail.Count == 0 || (Body.Position - Trail[^1]).Magnitude > 0.5)
            {
                Trail.Add(Body.Position);
                if (Trail.Count > MaxTrailLength)
                {
                    Trail.RemoveAt(0);
                }
            }

            // Check for decay
            if (!Properties.IsStable && !HasDecayed)
            {
                // Exponential decay: P(t) = e^(-t/τ) where τ = halflife/ln(2)
                double decayConstant = Math.Log(2) / Properties.HalfLife;
                double decayProbability = 1.0 - Math.Exp(-decayConstant * deltaTime);

                if (Rng.NextDouble() < decayProbability)
                {
                    Decay();
                }
            }
        }

        private void Decay()
        {
            HasDecayed = true;

            // Simple decay models (not fully realistic, but demonstrates the concept)
            switch (Properties.Type)
            {
                case ParticleType.Neutron:
                    // n → p + e⁻ + antineutrino (beta decay)
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Proton));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Electron));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Neutrino));
                    break;

                case ParticleType.Higgs:
                    // H → γ + γ (photon pair) or H → b + b̄
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                    break;

                case ParticleType.Quark:
                    // Quark → hadronization (simplified)
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Pion));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                    break;

                case ParticleType.Muon:
                    // μ⁻ → e⁻ + ν̄ₑ + νμ
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Electron));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Neutrino));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Neutrino));
                    break;

                case ParticleType.Pion:
                    // π⁰ → γ + γ
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                    DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                    break;
            }
        }

        private Particle2D CreateDecayProduct(ParticleType type)
        {
            // Random direction for decay products
            double angle = Rng.NextDouble() * Math.PI * 2;
            Vector2D direction = Vector2D.FromAngle(angle);

            // Conservation of momentum (simplified - lower velocity than parent)
            double speed = Body.Velocity.Magnitude * 0.5;
            Vector2D velocity = direction * speed;

            return new Particle2D(type, Body.Position, velocity);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double GetKineticEnergy()
        {
            return 0.5 * Properties.Mass * Body.Velocity.MagnitudeSquared;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double GetMomentum()
        {
            return Properties.Mass * Body.Velocity.Magnitude;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector2D GetMomentumVector()
        {
            return Body.Velocity * Properties.Mass;
        }

        /// <summary>
        /// Gets the relativistic gamma factor (Lorentz factor).
        /// </summary>
        public double GetLorentzFactor(double speedOfLight = 100.0)
        {
            double v = Body.Velocity.Magnitude;
            double beta = v / speedOfLight;
            if (beta >= 1.0) return double.PositiveInfinity;
            return 1.0 / Math.Sqrt(1.0 - beta * beta);
        }

        /// <summary>
        /// Gets the relativistic mass.
        /// </summary>
        public double GetRelativisticMass(double speedOfLight = 100.0)
        {
            return Properties.Mass * GetLorentzFactor(speedOfLight);
        }

        public override string ToString()
            => $"{Properties.Symbol} @ {Body.Position} v={Body.Velocity.Magnitude:F1}";
    }

    /// <summary>
    /// Particle beam - collection of particles moving together.
    /// </summary>
    public class ParticleBeam
    {
        public List<Particle2D> Particles { get; } = new();
        public string Name { get; set; } = "Beam";
        public bool IsClockwise { get; set; }

        public int Count => Particles.Count;

        public double TotalEnergy
        {
            get
            {
                double total = 0;
                foreach (var p in Particles)
                    total += p.GetKineticEnergy();
                return total;
            }
        }

        public double TotalMomentum
        {
            get
            {
                var total = Vector2D.Zero;
                foreach (var p in Particles)
                    total += p.GetMomentumVector();
                return total.Magnitude;
            }
        }

        public Vector2D AveragePosition
        {
            get
            {
                if (Particles.Count == 0) return Vector2D.Zero;
                var sum = Vector2D.Zero;
                foreach (var p in Particles)
                    sum += p.Body.Position;
                return sum / Particles.Count;
            }
        }

        public void Add(Particle2D particle)
        {
            Particles.Add(particle);
        }

        public void Clear()
        {
            Particles.Clear();
        }

        public void RemoveDecayed()
        {
            Particles.RemoveAll(p => p.HasDecayed);
        }
    }
}
