using Artemis.Physics2D;
using Artemis.Physics2D.Particles;

namespace ParticleColliderDemo;

/// <summary>
/// Console-friendly wrapper around Detector2D.
/// </summary>
public class Detector
{
    private readonly Detector2D _detector;

    public DetectorType2D Type => _detector.Type;
    public string Name => _detector.Name;
    public Vector2D Position => _detector.Position;
    public double Radius => _detector.Radius;
    public ConsoleColor Color => Type.GetConsoleColor();
    public RigidBody2D? TriggerZone => _detector.TriggerZone;

    public double TotalEnergyDetected => _detector.TotalEnergyDetected;
    public int ParticleCount => _detector.ParticleCount;

    public Detector(DetectorType2D type, string name, Vector2D position, double radius)
    {
        _detector = new Detector2D(type, name, position, radius);
    }

    public void RecordDetection(Particle particle, double simulationTime)
    {
        // Create a temporary Particle2D for detection (the wrapper's inner particle)
        // This is a bit of a workaround for the adapter pattern
        var p2d = new Particle2D(particle.Properties.Type, particle.Body.Position, particle.Body.Velocity);
        _detector.RecordDetection(p2d, simulationTime);
    }

    public void Clear()
    {
        _detector.Clear();
    }

    public string GetSummary()
    {
        return _detector.GetSummary();
    }
}

// Re-export enums for convenience
public enum DetectorType
{
    Tracker = DetectorType2D.Tracker,
    Calorimeter = DetectorType2D.Calorimeter,
    MuonChamber = DetectorType2D.MuonChamber,
    InteractionPoint = DetectorType2D.InteractionPoint
}
