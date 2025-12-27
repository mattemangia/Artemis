using ArtemisEngine;

namespace ParticleColliderDemo;

/// <summary>
/// Types of particles in the collider
/// Based on real particle physics
/// </summary>
public enum ParticleType
{
    Proton,      // Massive, positively charged
    Electron,    // Light, negatively charged
    Positron,    // Antimatter electron
    Neutron,     // Neutral, heavy
    Photon,      // Massless, light speed
    Higgs,       // Heavy boson (discovered at LHC!)
    Quark        // Fundamental particle
}

/// <summary>
/// Particle properties based on real physics
/// Masses in atomic mass units (amu), scaled for visualization
/// </summary>
public class ParticleProperties
{
    public ParticleType Type { get; set; }
    public float Mass { get; set; }          // Relative mass
    public int Charge { get; set; }           // Electric charge (-1, 0, +1)
    public ConsoleColor Color { get; set; }
    public required string Symbol { get; set; }
    public float Energy { get; set; }         // Kinetic energy (GeV scaled)
    public bool IsStable { get; set; }        // Does it decay?
    public float HalfLife { get; set; }       // Decay time (seconds)

    public static ParticleProperties Create(ParticleType type)
    {
        return type switch
        {
            ParticleType.Proton => new ParticleProperties
            {
                Type = ParticleType.Proton,
                Mass = 1.0f,
                Charge = +1,
                Color = ConsoleColor.Red,
                Symbol = "p⁺",
                IsStable = true,
                HalfLife = float.PositiveInfinity
            },
            ParticleType.Electron => new ParticleProperties
            {
                Type = ParticleType.Electron,
                Mass = 0.0005f, // ~1/2000 of proton
                Charge = -1,
                Color = ConsoleColor.Blue,
                Symbol = "e⁻",
                IsStable = true,
                HalfLife = float.PositiveInfinity
            },
            ParticleType.Positron => new ParticleProperties
            {
                Type = ParticleType.Positron,
                Mass = 0.0005f,
                Charge = +1,
                Color = ConsoleColor.Cyan,
                Symbol = "e⁺",
                IsStable = true,
                HalfLife = float.PositiveInfinity
            },
            ParticleType.Neutron => new ParticleProperties
            {
                Type = ParticleType.Neutron,
                Mass = 1.0f,
                Charge = 0,
                Color = ConsoleColor.Gray,
                Symbol = "n⁰",
                IsStable = false,
                HalfLife = 881.5f // seconds
            },
            ParticleType.Photon => new ParticleProperties
            {
                Type = ParticleType.Photon,
                Mass = 0.0001f, // Nearly massless
                Charge = 0,
                Color = ConsoleColor.Yellow,
                Symbol = "γ",
                IsStable = true,
                HalfLife = float.PositiveInfinity
            },
            ParticleType.Higgs => new ParticleProperties
            {
                Type = ParticleType.Higgs,
                Mass = 125.0f, // GeV/c² (scaled)
                Charge = 0,
                Color = ConsoleColor.Magenta,
                Symbol = "H⁰",
                IsStable = false,
                HalfLife = 1.6e-22f // Extremely short!
            },
            ParticleType.Quark => new ParticleProperties
            {
                Type = ParticleType.Quark,
                Mass = 0.3f,
                Charge = -1,
                Color = ConsoleColor.Green,
                Symbol = "q",
                IsStable = false,
                HalfLife = 1e-23f
            },
            _ => throw new ArgumentException($"Unknown particle type: {type}")
        };
    }
}

/// <summary>
/// Wrapper around RigidBody with particle-specific properties
/// </summary>
public class Particle
{
    public RigidBody Body { get; set; }
    public ParticleProperties Properties { get; set; }
    public float Age { get; set; } = 0;
    public bool HasDecayed { get; set; } = false;
    public List<Particle> DecayProducts { get; set; } = new();
    public Vector2 InitialPosition { get; set; }
    public float TotalEnergyDeposited { get; set; } = 0; // Energy in detectors

    // Trail for visualization
    public List<Vector2> Trail { get; set; } = new();
    public const int MaxTrailLength = 50;

    public Particle(ParticleType type, Vector2 position, Vector2 velocity)
    {
        Properties = ParticleProperties.Create(type);

        // Create physics body (very small particles)
        float visualRadius = MathF.Max(0.1f, Properties.Mass * 0.15f);
        Body = new RigidBody(position, Properties.Mass, new CircleShape(visualRadius));
        Body.Velocity = velocity;
        Body.Restitution = 1.0f; // Perfectly elastic (ideal)
        Body.Friction = 0.0f;     // No friction in vacuum
        Body.UseCCD = true;       // Essential for high-speed particles!
        Body.UserData = this;     // Link back to particle

        InitialPosition = position;

        // Calculate initial kinetic energy
        Properties.Energy = 0.5f * Properties.Mass * velocity.LengthSquared;
    }

    public void Update(float deltaTime)
    {
        Age += deltaTime;

        // Record trail
        if (Trail.Count == 0 || (Body.Position - Trail[^1]).Length > 0.5f)
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
            float decayConstant = MathF.Log(2) / Properties.HalfLife;
            float decayProbability = 1.0f - MathF.Exp(-decayConstant * deltaTime);

            if (Random.Shared.NextSingle() < decayProbability)
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
                break;

            case ParticleType.Higgs:
                // H → γ + γ (photon pair)
                DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                break;

            case ParticleType.Quark:
                // Quark → hadronization (simplified)
                DecayProducts.Add(CreateDecayProduct(ParticleType.Photon));
                break;
        }
    }

    private Particle CreateDecayProduct(ParticleType type)
    {
        // Random direction for decay products
        float angle = Random.Shared.NextSingle() * MathF.PI * 2;
        Vector2 direction = new Vector2(MathF.Cos(angle), MathF.Sin(angle));

        // Lower velocity than parent (energy conservation, but simplified)
        Vector2 velocity = direction * Body.Velocity.Length * 0.5f;

        return new Particle(type, Body.Position, velocity);
    }

    public float GetKineticEnergy()
    {
        return 0.5f * Properties.Mass * Body.Velocity.LengthSquared;
    }

    public float GetMomentum()
    {
        return Properties.Mass * Body.Velocity.Length;
    }
}
