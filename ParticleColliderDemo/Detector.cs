using ArtemisEngine;

namespace ParticleColliderDemo;

/// <summary>
/// Types of detectors in a particle collider (based on LHC)
/// </summary>
public enum DetectorType
{
    Tracker,           // Tracks particle paths
    Calorimeter,       // Measures energy
    MuonChamber,       // Detects muons
    InteractionPoint   // Where collisions happen
}

/// <summary>
/// Detector that records particles passing through
/// Similar to ATLAS or CMS detectors at LHC
/// </summary>
public class Detector
{
    public DetectorType Type { get; set; }
    public string Name { get; set; }
    public Vector2 Position { get; set; }
    public float Radius { get; set; }
    public ConsoleColor Color { get; set; }
    public RigidBody TriggerZone { get; set; }

    // Detection data
    public List<ParticleDetection> Detections { get; private set; } = new();
    public float TotalEnergyDetected { get; private set; } = 0;
    public int ParticleCount { get; private set; } = 0;

    public Detector(DetectorType type, string name, Vector2 position, float radius)
    {
        Type = type;
        Name = name;
        Position = position;
        Radius = radius;

        Color = type switch
        {
            DetectorType.Tracker => ConsoleColor.Cyan,
            DetectorType.Calorimeter => ConsoleColor.Yellow,
            DetectorType.MuonChamber => ConsoleColor.Magenta,
            DetectorType.InteractionPoint => ConsoleColor.Red,
            _ => ConsoleColor.White
        };

        // Create trigger zone (no physics, just detection)
        TriggerZone = new RigidBody(position, 0, new CircleShape(radius));
        TriggerZone.IsStatic = true;
        TriggerZone.IsTrigger = true;
        TriggerZone.UserData = this;
    }

    public void RecordDetection(Particle particle, float simulationTime)
    {
        var detection = new ParticleDetection
        {
            Time = simulationTime,
            ParticleType = particle.Properties.Type,
            Position = particle.Body.Position,
            Velocity = particle.Body.Velocity,
            Energy = particle.GetKineticEnergy(),
            Momentum = particle.GetMomentum(),
            Charge = particle.Properties.Charge
        };

        Detections.Add(detection);
        TotalEnergyDetected += detection.Energy;
        ParticleCount++;
    }

    public void Clear()
    {
        Detections.Clear();
        TotalEnergyDetected = 0;
        ParticleCount = 0;
    }

    public string GetSummary()
    {
        if (ParticleCount == 0)
            return $"{Name}: No detections";

        return $"{Name}: {ParticleCount} particles, {TotalEnergyDetected:F1} GeV";
    }
}

/// <summary>
/// Record of a particle detection
/// </summary>
public class ParticleDetection
{
    public float Time { get; set; }
    public ParticleType ParticleType { get; set; }
    public Vector2 Position { get; set; }
    public Vector2 Velocity { get; set; }
    public float Energy { get; set; }
    public float Momentum { get; set; }
    public int Charge { get; set; }
}

/// <summary>
/// Collision event between particles
/// Represents a high-energy collision like those at the LHC
/// </summary>
public class ParticleCollisionEvent
{
    public float Time { get; set; }
    public Particle Particle1 { get; set; }
    public Particle Particle2 { get; set; }
    public Vector2 CollisionPoint { get; set; }
    public float TotalEnergyBefore { get; set; }
    public float TotalEnergyAfter { get; set; }
    public float CenterOfMassEnergy { get; set; }
    public List<Particle> Products { get; set; } = new();

    public ParticleCollisionEvent(Particle p1, Particle p2, Vector2 point, float time)
    {
        Time = time;
        Particle1 = p1;
        Particle2 = p2;
        CollisionPoint = point;

        // Calculate center-of-mass energy (like √s at LHC)
        float e1 = p1.GetKineticEnergy();
        float e2 = p2.GetKineticEnergy();
        TotalEnergyBefore = e1 + e2;

        // Simplified COM energy calculation
        CenterOfMassEnergy = MathF.Sqrt(2 * p1.Properties.Mass * p2.Properties.Mass +
                                        2 * e1 * e2);
    }

    public void AnalyzeProducts(List<Particle> allParticles)
    {
        // Find particles created near collision point recently
        foreach (var p in allParticles)
        {
            if ((p.Body.Position - CollisionPoint).Length < 2.0f &&
                p.Age < 0.1f &&
                p != Particle1 && p != Particle2)
            {
                Products.Add(p);
            }
        }

        TotalEnergyAfter = Products.Sum(p => p.GetKineticEnergy());
    }

    public string GetSummary()
    {
        string p1Type = Particle1.Properties.Symbol;
        string p2Type = Particle2.Properties.Symbol;

        return $"{p1Type} + {p2Type} → {Products.Count} products, " +
               $"√s = {CenterOfMassEnergy:F1} GeV, " +
               $"E_in = {TotalEnergyBefore:F1} GeV";
    }
}
