using Artemis.Physics2D;
using Artemis.Physics2D.Particles;

namespace ParticleColliderDemo;

/// <summary>
/// Console-friendly wrapper around Particle2D with ConsoleColor support.
/// </summary>
public class Particle
{
    private readonly Particle2D _particle;

    public RigidBody2D Body => _particle.Body;
    public ParticleProperties Properties => _particle.Properties;
    public double Age => _particle.Age;
    public bool HasDecayed => _particle.HasDecayed;
    public List<Particle2D> DecayProducts => _particle.DecayProducts;
    public Vector2D InitialPosition => _particle.InitialPosition;
    public double TotalEnergyDeposited
    {
        get => _particle.TotalEnergyDeposited;
        set => _particle.TotalEnergyDeposited = value;
    }

    public List<Vector2D> Trail => _particle.Trail;
    public int MaxTrailLength
    {
        get => _particle.MaxTrailLength;
        set => _particle.MaxTrailLength = value;
    }

    public Particle(ParticleType type, Vector2D position, Vector2D velocity)
    {
        _particle = new Particle2D(type, position, velocity);
        // Override UserData to point to this wrapper instead of the underlying Particle2D
        _particle.Body.UserData = this;
    }

    public void Update(double deltaTime)
    {
        _particle.Update(deltaTime);
    }

    public double GetKineticEnergy() => _particle.GetKineticEnergy();
    public double GetMomentum() => _particle.GetMomentum();
    public Vector2D GetMomentumVector() => _particle.GetMomentumVector();
}

/// <summary>
/// Extension methods for console color support.
/// </summary>
public static class ParticleExtensions
{
    public static ConsoleColor GetConsoleColor(this ParticleProperties props)
    {
        return props.Type switch
        {
            ParticleType.Proton => ConsoleColor.Red,
            ParticleType.Electron => ConsoleColor.Blue,
            ParticleType.Positron => ConsoleColor.Cyan,
            ParticleType.Neutron => ConsoleColor.Gray,
            ParticleType.Photon => ConsoleColor.Yellow,
            ParticleType.Higgs => ConsoleColor.Magenta,
            ParticleType.Quark => ConsoleColor.Green,
            ParticleType.Muon => ConsoleColor.DarkYellow,
            ParticleType.Neutrino => ConsoleColor.White,
            ParticleType.Pion => ConsoleColor.DarkMagenta,
            _ => ConsoleColor.White
        };
    }

    public static ConsoleColor GetConsoleColor(this DetectorType2D type)
    {
        return type switch
        {
            DetectorType2D.Tracker => ConsoleColor.Cyan,
            DetectorType2D.Calorimeter => ConsoleColor.Yellow,
            DetectorType2D.MuonChamber => ConsoleColor.Magenta,
            DetectorType2D.InteractionPoint => ConsoleColor.Red,
            DetectorType2D.TriggerSystem => ConsoleColor.Green,
            DetectorType2D.VertexDetector => ConsoleColor.DarkYellow,
            _ => ConsoleColor.White
        };
    }
}
