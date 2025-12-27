using ArtemisEngine;

namespace ParticleColliderDemo;

/// <summary>
/// The particle accelerator ring (like the LHC's 27km ring)
/// Uses magnetic field effectors to keep particles in orbit
/// </summary>
public class Accelerator
{
    public Vector2 Center { get; set; }
    public float Radius { get; set; }
    public float BeamEnergy { get; set; } // Energy in GeV

    // Two counter-rotating beams (like LHC)
    public List<Particle> Beam1 { get; private set; } = new(); // Clockwise
    public List<Particle> Beam2 { get; private set; } = new(); // Counter-clockwise

    // Magnetic confinement
    private RadialForceEffector _magneticField;

    public Accelerator(Vector2 center, float radius, float beamEnergy)
    {
        Center = center;
        Radius = radius;
        BeamEnergy = beamEnergy;

        // Create magnetic field to keep particles in orbit
        _magneticField = new RadialForceEffector(
            center: center,
            radius: radius + 2,
            strength: -50f, // Attractive force toward center
            falloff: RadialFalloff.Linear
        );
    }

    public RadialForceEffector GetMagneticField() => _magneticField;

    /// <summary>
    /// Inject particles into beam 1 (clockwise)
    /// </summary>
    public void InjectBeam1(ParticleType type, int count)
    {
        for (int i = 0; i < count; i++)
        {
            float angle = (float)i / count * MathF.PI * 2;
            Vector2 position = Center + new Vector2(
                MathF.Cos(angle) * Radius,
                MathF.Sin(angle) * Radius
            );

            // Velocity tangent to circle (clockwise)
            Vector2 tangent = new Vector2(-MathF.Sin(angle), MathF.Cos(angle));
            float speed = CalculateBeamSpeed(type);
            Vector2 velocity = tangent * speed;

            var particle = new Particle(type, position, velocity);
            Beam1.Add(particle);
        }
    }

    /// <summary>
    /// Inject particles into beam 2 (counter-clockwise)
    /// </summary>
    public void InjectBeam2(ParticleType type, int count)
    {
        for (int i = 0; i < count; i++)
        {
            float angle = (float)i / count * MathF.PI * 2;
            Vector2 position = Center + new Vector2(
                MathF.Cos(angle) * Radius,
                MathF.Sin(angle) * Radius
            );

            // Velocity tangent to circle (counter-clockwise)
            Vector2 tangent = new Vector2(MathF.Sin(angle), -MathF.Cos(angle));
            float speed = CalculateBeamSpeed(type);
            Vector2 velocity = tangent * speed;

            var particle = new Particle(type, position, velocity);
            Beam2.Add(particle);
        }
    }

    /// <summary>
    /// Calculate particle speed based on beam energy
    /// E = 0.5 * m * v^2, so v = sqrt(2*E/m)
    /// </summary>
    private float CalculateBeamSpeed(ParticleType type)
    {
        var props = ParticleProperties.Create(type);
        // Avoid division by zero for massless particles
        float mass = MathF.Max(props.Mass, 0.0001f);
        return MathF.Sqrt(2 * BeamEnergy / mass);
    }

    /// <summary>
    /// Get interaction points where beams should cross
    /// At LHC, there are 4 main interaction points (ATLAS, CMS, ALICE, LHCb)
    /// </summary>
    public List<Vector2> GetInteractionPoints()
    {
        return new List<Vector2>
        {
            Center + new Vector2(Radius, 0),      // IP1 (ATLAS-like)
            Center + new Vector2(0, Radius),      // IP2
            Center + new Vector2(-Radius, 0),     // IP3 (CMS-like)
            Center + new Vector2(0, -Radius)      // IP4
        };
    }

    /// <summary>
    /// Apply focusing to keep particles in beam
    /// (simplified model of quadrupole magnets)
    /// </summary>
    public void ApplyBeamFocusing(Particle particle, float deltaTime)
    {
        // Calculate distance from ideal orbit
        float distanceFromCenter = (particle.Body.Position - Center).Length;
        float deviation = distanceFromCenter - Radius;

        if (MathF.Abs(deviation) > 0.1f)
        {
            // Apply restoring force toward ideal orbit
            Vector2 toCenter = (Center - particle.Body.Position).Normalized;
            Vector2 restoringForce = toCenter * (deviation * 10f);
            particle.Body.ApplyForce(restoringForce * deltaTime);
        }
    }

    /// <summary>
    /// Get all particles in the accelerator
    /// </summary>
    public IEnumerable<Particle> GetAllParticles()
    {
        return Beam1.Concat(Beam2);
    }

    /// <summary>
    /// Clear all beams
    /// </summary>
    public void ClearBeams()
    {
        Beam1.Clear();
        Beam2.Clear();
    }

    /// <summary>
    /// Get beam statistics
    /// </summary>
    public (int beam1Count, int beam2Count, float avgEnergy) GetStatistics()
    {
        int b1 = Beam1.Count;
        int b2 = Beam2.Count;
        float avgE = 0;

        var all = GetAllParticles().ToList();
        if (all.Count > 0)
        {
            avgE = all.Average(p => p.GetKineticEnergy());
        }

        return (b1, b2, avgE);
    }
}

/// <summary>
/// Experiment coordinator - manages the collision experiment
/// </summary>
public class CollisionExperiment
{
    public string Name { get; set; }
    public Accelerator Accelerator { get; set; }
    public List<Detector> Detectors { get; private set; } = new();
    public List<ParticleCollisionEvent> CollisionEvents { get; private set; } = new();

    public float SimulationTime { get; private set; } = 0;
    public int TotalCollisions { get; private set; } = 0;

    // Data export
    private SimulationDataExporter _dataExporter = new();
    private EnergyTracker _energyTracker = new();

    public CollisionExperiment(string name, Accelerator accelerator)
    {
        Name = name;
        Accelerator = accelerator;
        SetupDetectors();
    }

    private void SetupDetectors()
    {
        var interactionPoints = Accelerator.GetInteractionPoints();

        // Create detector array at each interaction point
        for (int i = 0; i < interactionPoints.Count; i++)
        {
            var ip = interactionPoints[i];

            // Inner tracker
            Detectors.Add(new Detector(
                DetectorType.Tracker,
                $"Tracker_{i + 1}",
                ip,
                3.0f
            ));

            // Calorimeter (outer layer)
            Detectors.Add(new Detector(
                DetectorType.Calorimeter,
                $"Calorimeter_{i + 1}",
                ip,
                5.0f
            ));

            // Interaction point marker
            Detectors.Add(new Detector(
                DetectorType.InteractionPoint,
                $"IP{i + 1}",
                ip,
                0.5f
            ));
        }
    }

    public void Update(float deltaTime, PhysicsWorld world)
    {
        SimulationTime += deltaTime;

        // Update energy tracking
        var allBodies = world.Bodies.Where(b => b.UserData is Particle).ToList();
        _energyTracker.UpdateEnergy(allBodies);
    }

    public void RecordCollision(Particle p1, Particle p2, Vector2 point)
    {
        var collisionEvent = new ParticleCollisionEvent(p1, p2, point, SimulationTime);
        CollisionEvents.Add(collisionEvent);
        TotalCollisions++;
    }

    public string GetEnergyReport()
    {
        return $"KE: {_energyTracker.KineticEnergy:F1} GeV | " +
               $"Total: {_energyTracker.TotalEnergy:F1} GeV | " +
               $"Rot: {_energyTracker.RotationalEnergy:F1} GeV";
    }

    public void ExportData(string filename)
    {
        string csv = _dataExporter.ExportToCSV();
        File.WriteAllText(filename, csv);
    }
}
