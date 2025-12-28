using Artemis.Physics2D;
using Artemis.Physics2D.Particles;
using ParticleEnergyTracker2D = Artemis.Physics2D.Particles.EnergyTracker2D;
using ParticleSimulationDataExporter2D = Artemis.Physics2D.Particles.SimulationDataExporter2D;

namespace ParticleColliderDemo;

/// <summary>
/// Console-friendly wrapper around Accelerator2D.
/// </summary>
public class Accelerator
{
    private readonly Accelerator2D _accelerator;

    public Vector2D Center => _accelerator.Center;
    public double Radius => _accelerator.Radius;
    public double BeamEnergy => _accelerator.BeamEnergy;

    public List<Particle> Beam1 { get; } = new();
    public List<Particle> Beam2 { get; } = new();

    public Accelerator(Vector2D center, double radius, double beamEnergy)
    {
        _accelerator = new Accelerator2D(center, radius, beamEnergy);
    }

    public RadialForceEffector2D GetMagneticField() => _accelerator.GetMagneticField();

    public void InjectBeam1(ParticleType type, int count)
    {
        for (int i = 0; i < count; i++)
        {
            double angle = (double)i / count * Math.PI * 2;
            Vector2D position = Center + Vector2D.FromAngle(angle, Radius);

            // Velocity tangent to circle (clockwise)
            Vector2D tangent = Vector2D.FromAngle(angle - Math.PI / 2);
            double speed = CalculateBeamSpeed(type);
            Vector2D velocity = tangent * speed;

            var particle = new Particle(type, position, velocity);
            Beam1.Add(particle);
        }
    }

    public void InjectBeam2(ParticleType type, int count)
    {
        for (int i = 0; i < count; i++)
        {
            double angle = (double)i / count * Math.PI * 2;
            Vector2D position = Center + Vector2D.FromAngle(angle, Radius);

            // Velocity tangent to circle (counter-clockwise)
            Vector2D tangent = Vector2D.FromAngle(angle + Math.PI / 2);
            double speed = CalculateBeamSpeed(type);
            Vector2D velocity = tangent * speed;

            var particle = new Particle(type, position, velocity);
            Beam2.Add(particle);
        }
    }

    private double CalculateBeamSpeed(ParticleType type)
    {
        var props = ParticleProperties.Create(type);
        double mass = Math.Max(props.Mass, 0.0001);
        return Math.Sqrt(2 * BeamEnergy / mass);
    }

    public List<Vector2D> GetInteractionPoints()
    {
        return new List<Vector2D>
        {
            Center + new Vector2D(Radius, 0),
            Center + new Vector2D(0, Radius),
            Center + new Vector2D(-Radius, 0),
            Center + new Vector2D(0, -Radius)
        };
    }

    public void ApplyBeamFocusing(Particle particle, double deltaTime)
    {
        double distanceFromCenter = (particle.Body.Position - Center).Magnitude;
        double deviation = distanceFromCenter - Radius;

        if (Math.Abs(deviation) > 0.1)
        {
            Vector2D toCenter = (Center - particle.Body.Position).Normalized;
            Vector2D restoringForce = toCenter * (deviation * 10.0);
            particle.Body.ApplyForce(restoringForce * deltaTime);
        }
    }

    public void Update(double deltaTime)
    {
        // 1. Apply RF Acceleration
        _accelerator.ApplyRFAccelerationParallel(deltaTime);

        // 2. Adjust Magnetic Field Strength (Synchrotron principle)
        // B = p / (q * r)  where p is momentum
        // Force F = qvB = mv^2/r -> B = mv/qr = p/qr

        // We calculate required B for the average energy/momentum of the beam
        var (b1, b2, avgE) = GetStatistics();
        if (avgE > 0)
        {
            // Approximate average speed/momentum
            // For simplicity, let's just scale B with Energy roughly or use a strong enough force
            // F_centripetal = m * v^2 / r
            // Kinetic Energy K = 0.5 * m * v^2 -> v^2 = 2K/m
            // F = m * (2K/m) / r = 2K / r

            // So we need a force of magnitude 2K/r pointing to center.
            // The RadialForceEffector applies a force.
            // We set its strength to -2K/r (negative for attraction)

            double requiredForce = -(2.0 * avgE) / Radius;

            // Add a safety factor to keep them tight
            requiredForce *= 1.2;

            // Limit to max reasonable force to prevent numerical explosion
             if (Math.Abs(requiredForce) > 10000) requiredForce = -10000;
             if (Math.Abs(requiredForce) < 50) requiredForce = -50;

            _accelerator.SetMagneticFieldStrength(requiredForce);
        }
    }

    public IEnumerable<Particle> GetAllParticles()
    {
        return Beam1.Concat(Beam2);
    }

    public void ClearBeams()
    {
        Beam1.Clear();
        Beam2.Clear();
    }

    public (int beam1Count, int beam2Count, double avgEnergy) GetStatistics()
    {
        int b1 = Beam1.Count;
        int b2 = Beam2.Count;
        double avgE = 0;

        var all = GetAllParticles().ToList();
        if (all.Count > 0)
        {
            avgE = all.Average(p => p.GetKineticEnergy());
        }

        return (b1, b2, avgE);
    }
}

/// <summary>
/// Console-friendly experiment coordinator.
/// </summary>
public class CollisionExperiment
{
    public string Name { get; set; }
    public Accelerator Accelerator { get; }
    public List<Detector> Detectors { get; } = new();
    public List<ParticleCollisionEvent> CollisionEvents { get; } = new();

    public double SimulationTime { get; private set; }
    public int TotalCollisions { get; private set; }

    private ParticleEnergyTracker2D _energyTracker = new();
    private ParticleSimulationDataExporter2D _dataExporter = new();

    public CollisionExperiment(string name, Accelerator accelerator)
    {
        Name = name;
        Accelerator = accelerator;
        SetupDetectors();
    }

    private void SetupDetectors()
    {
        var interactionPoints = Accelerator.GetInteractionPoints();

        for (int i = 0; i < interactionPoints.Count; i++)
        {
            var ip = interactionPoints[i];

            Detectors.Add(new Detector(
                DetectorType2D.Tracker,
                $"Tracker_{i + 1}",
                ip,
                3.0
            ));

            Detectors.Add(new Detector(
                DetectorType2D.Calorimeter,
                $"Calorimeter_{i + 1}",
                ip,
                5.0
            ));

            Detectors.Add(new Detector(
                DetectorType2D.InteractionPoint,
                $"IP{i + 1}",
                ip,
                0.5
            ));
        }
    }

    public void Update(double deltaTime, PhysicsWorld2D world)
    {
        SimulationTime += deltaTime;

        var allBodies = world.Bodies
            .Where(b => b.UserData is Particle)
            .ToList();
        _energyTracker.UpdateEnergy(allBodies);
    }

    public void RecordCollision(Particle p1, Particle p2, Vector2D point)
    {
        var collisionEvent = new ParticleCollisionEvent(p1, p2, point, SimulationTime);
        CollisionEvents.Add(collisionEvent);
        TotalCollisions++;
    }

    public string GetEnergyReport()
    {
        return $"KE: {_energyTracker.KineticEnergy:F1} GeV | " +
               $"Total: {_energyTracker.TotalEnergy:F1} GeV";
    }

    public void ExportData(string filename)
    {
        var csv = new System.Text.StringBuilder();
        csv.AppendLine("Time,ParticleType,Energy,Collision");

        foreach (var evt in CollisionEvents)
        {
            csv.AppendLine($"{evt.Time:F4},{evt.Particle1.Properties.Type},{evt.TotalEnergyBefore:F2},Yes");
        }

        File.WriteAllText(filename, csv.ToString());
    }
}

/// <summary>
/// Collision event between particles.
/// </summary>
public class ParticleCollisionEvent
{
    public double Time { get; set; }
    public Particle Particle1 { get; set; }
    public Particle Particle2 { get; set; }
    public Vector2D CollisionPoint { get; set; }
    public double TotalEnergyBefore { get; set; }
    public double TotalEnergyAfter { get; set; }
    public double CenterOfMassEnergy { get; set; }
    public List<Particle> Products { get; set; } = new();

    public ParticleCollisionEvent(Particle p1, Particle p2, Vector2D point, double time)
    {
        Time = time;
        Particle1 = p1;
        Particle2 = p2;
        CollisionPoint = point;

        double e1 = p1.GetKineticEnergy();
        double e2 = p2.GetKineticEnergy();
        TotalEnergyBefore = e1 + e2;

        CenterOfMassEnergy = Math.Sqrt(2 * p1.Properties.Mass * p2.Properties.Mass +
                                        2 * e1 * e2);
    }

    public void AnalyzeProducts(List<Particle> allParticles)
    {
        foreach (var p in allParticles)
        {
            if ((p.Body.Position - CollisionPoint).Magnitude < 2.0 &&
                p.Age < 0.1 &&
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
