using Artemis.Graphics;
using ArtemisEngine;
using OpenTK.Mathematics;
using OpenTK.Windowing.GraphicsLibraryFramework;
using ParticleColliderDemo;

namespace ParticleColliderDemo.OpenTK;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== ARTEMIS PARTICLE COLLIDER - OpenTK Edition ===");
        Console.WriteLine("Starting graphical simulation...");

        using var window = new ParticleColliderWindow(1280, 720, "Artemis Particle Collider - LHC Simulator");
        window.Run();
    }
}

public class ParticleColliderWindow : GraphicsWindow
{
    private Accelerator _accelerator = null!;
    private CollisionExperiment _experiment = null!;
    private List<(ParticleCollisionEvent evt, float timestamp)> _recentCollisions = new();

    // Colors for particle types
    private static readonly Dictionary<ParticleType, Color4> ParticleColors = new()
    {
        { ParticleType.Proton, new Color4(1f, 0.2f, 0.2f, 1f) },
        { ParticleType.Electron, new Color4(0.2f, 0.4f, 1f, 1f) },
        { ParticleType.Positron, new Color4(0.2f, 0.9f, 0.9f, 1f) },
        { ParticleType.Neutron, new Color4(0.6f, 0.6f, 0.6f, 1f) },
        { ParticleType.Photon, new Color4(1f, 1f, 0.2f, 1f) },
        { ParticleType.Higgs, new Color4(0.9f, 0.2f, 0.9f, 1f) },
        { ParticleType.Quark, new Color4(0.2f, 0.9f, 0.2f, 1f) },
    };

    public ParticleColliderWindow(int width, int height, string title)
        : base(width, height, title)
    {
        // Slow down for visibility
        TimeScale = 0.3f;
    }

    protected override void Initialize()
    {
        // No gravity in the accelerator (vacuum)
        World = new PhysicsWorld(new ArtemisEngine.Vector2(0, 0));
        World.UseAdvancedSolver = true;

        // Create accelerator with smaller energy for slower particles
        _accelerator = new Accelerator(
            center: new ArtemisEngine.Vector2(0, 0),
            radius: 15f,
            beamEnergy: 3f // Lower energy = slower particles
        );

        // Create experiment
        _experiment = new CollisionExperiment("ARTEMIS Particle Physics Experiment", _accelerator);

        // Add magnetic field
        World.AddAreaEffector(_accelerator.GetMagneticField());

        // Add detector trigger zones
        foreach (var detector in _experiment.Detectors)
        {
            World.AddBody(detector.TriggerZone);
        }

        // Setup collision events
        World.OnTriggerEnter += (sender, e) =>
        {
            if (e.BodyA.UserData is Detector detector && e.BodyB.UserData is Particle particle)
            {
                detector.RecordDetection(particle, _experiment.SimulationTime);
            }
            else if (e.BodyB.UserData is Detector detector2 && e.BodyA.UserData is Particle particle2)
            {
                detector2.RecordDetection(particle2, _experiment.SimulationTime);
            }
        };

        World.OnCollisionEnter += (sender, e) =>
        {
            if (e.BodyA.UserData is Particle p1 && e.BodyB.UserData is Particle p2)
            {
                ArtemisEngine.Vector2 collisionPoint = (p1.Body.Position + p2.Body.Position) * 0.5f;
                _experiment.RecordCollision(p1, p2, collisionPoint);

                var evt = new ParticleCollisionEvent(p1, p2, collisionPoint, _experiment.SimulationTime);
                _recentCollisions.Add((evt, _experiment.SimulationTime));

                // Create collision products
                if (Random.Shared.NextSingle() < 0.3f)
                {
                    CreateCollisionProducts(collisionPoint, p1, p2);
                }
            }
        };

        // Inject initial beams
        InjectInitialBeams();

        // Set camera
        CameraZoom = 25;
    }

    private void InjectInitialBeams()
    {
        _accelerator.InjectBeam1(ParticleType.Proton, 6);
        _accelerator.InjectBeam2(ParticleType.Proton, 6);

        foreach (var particle in _accelerator.GetAllParticles())
        {
            World.AddBody(particle.Body);
        }
    }

    private void InjectBeams()
    {
        ParticleType[] types = { ParticleType.Proton, ParticleType.Electron, ParticleType.Positron };
        var type1 = types[Random.Shared.Next(types.Length)];
        var type2 = types[Random.Shared.Next(types.Length)];

        _accelerator.InjectBeam1(type1, 4);
        _accelerator.InjectBeam2(type2, 4);

        foreach (var particle in _accelerator.Beam1.TakeLast(4).Concat(_accelerator.Beam2.TakeLast(4)))
        {
            World.AddBody(particle.Body);
        }
    }

    private void InjectSpecificBeam(ParticleType type1, ParticleType type2)
    {
        _accelerator.InjectBeam1(type1, 5);
        _accelerator.InjectBeam2(type2, 5);

        foreach (var particle in _accelerator.Beam1.TakeLast(5).Concat(_accelerator.Beam2.TakeLast(5)))
        {
            World.AddBody(particle.Body);
        }
    }

    private void CreateCollisionProducts(ArtemisEngine.Vector2 point, Particle p1, Particle p2)
    {
        float availableEnergy = p1.GetKineticEnergy() + p2.GetKineticEnergy();
        int productCount = Random.Shared.Next(2, 4);

        for (int i = 0; i < productCount; i++)
        {
            float angle = Random.Shared.NextSingle() * MathF.PI * 2;
            ArtemisEngine.Vector2 direction = new ArtemisEngine.Vector2(MathF.Cos(angle), MathF.Sin(angle));

            float energyFraction = Random.Shared.NextSingle() * 0.4f;
            float speed = MathF.Sqrt(2 * availableEnergy * energyFraction);

            ParticleType productType = Random.Shared.NextSingle() < 0.5f
                ? ParticleType.Photon
                : ParticleType.Electron;

            var product = new Particle(productType, point, direction * speed);
            World.AddBody(product.Body);
        }
    }

    private void ClearExperiment()
    {
        var particleBodies = World.Bodies
            .Where(b => b.UserData is Particle)
            .ToList();

        foreach (var body in particleBodies)
        {
            World.RemoveBody(body);
        }

        _accelerator.ClearBeams();

        foreach (var detector in _experiment.Detectors)
        {
            detector.Clear();
        }

        _recentCollisions.Clear();
    }

    protected override void HandleInput()
    {
        // Inject beams
        if (KeyboardState.IsKeyPressed(Keys.Space))
        {
            InjectBeams();
        }

        // Clear experiment
        if (KeyboardState.IsKeyPressed(Keys.C))
        {
            ClearExperiment();
        }

        // Specific beam types
        if (KeyboardState.IsKeyPressed(Keys.D1))
        {
            InjectSpecificBeam(ParticleType.Proton, ParticleType.Proton);
        }
        if (KeyboardState.IsKeyPressed(Keys.D2))
        {
            InjectSpecificBeam(ParticleType.Electron, ParticleType.Positron);
        }
        if (KeyboardState.IsKeyPressed(Keys.D3))
        {
            InjectSpecificBeam(ParticleType.Proton, ParticleType.Electron);
        }

        // Export data
        if (KeyboardState.IsKeyPressed(Keys.E))
        {
            try
            {
                string filename = $"particle_collision_data_{DateTime.Now:yyyyMMdd_HHmmss}.csv";
                _experiment.ExportData(filename);
                Console.WriteLine($"Data exported to: {filename}");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Export failed: {ex.Message}");
            }
        }
    }

    protected override void UpdatePhysics(float deltaTime)
    {
        // Update particles
        var allParticles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .ToList();

        foreach (var particle in allParticles)
        {
            particle.Update(deltaTime);

            // Apply beam focusing
            _accelerator.ApplyBeamFocusing(particle, deltaTime);

            // Remove decayed particles
            if (particle.HasDecayed)
            {
                World.RemoveBody(particle.Body);

                foreach (var product in particle.DecayProducts)
                {
                    World.AddBody(product.Body);
                }
            }
        }

        // Update experiment time
        _experiment.Update(deltaTime, World);

        // Clean up old collision events
        _recentCollisions.RemoveAll(c => _experiment.SimulationTime - c.timestamp > 2.0f);
    }

    protected override void Render()
    {
        // Draw accelerator ring
        DrawRing(new ArtemisEngine.Vector2(0, 0), _accelerator.Radius, new Color4(0.3f, 0.3f, 0.4f, 1f), 3f);

        // Draw interaction points
        foreach (var ip in _accelerator.GetInteractionPoints())
        {
            DrawCircle(ip, 0.8f, new Color4(0.8f, 0.2f, 0.2f, 1f), true);
            DrawCircle(ip, 1.2f, new Color4(1f, 0.3f, 0.3f, 0.5f), false);
        }

        // Draw detectors
        foreach (var detector in _experiment.Detectors)
        {
            Color4 detectorColor = detector.Type switch
            {
                DetectorType.Tracker => new Color4(0.2f, 0.5f, 0.8f, 0.3f),
                DetectorType.Calorimeter => new Color4(0.5f, 0.8f, 0.2f, 0.3f),
                DetectorType.MuonChamber => new Color4(0.8f, 0.5f, 0.2f, 0.3f),
                DetectorType.InteractionPoint => new Color4(1f, 0.2f, 0.2f, 0.5f),
                _ => new Color4(0.5f, 0.5f, 0.5f, 0.3f)
            };

            DrawCircle(detector.Position, detector.Radius, detectorColor, false);
        }

        // Get all particles
        var allParticles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .ToList();

        // Draw particle trails
        foreach (var particle in allParticles)
        {
            if (particle.Trail.Count > 1)
            {
                Color4 particleColor = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);
                Color4 trailStart = new Color4(particleColor.R * 0.3f, particleColor.G * 0.3f, particleColor.B * 0.3f, 0.2f);
                Color4 trailEnd = new Color4(particleColor.R, particleColor.G, particleColor.B, 0.8f);

                DrawTrail(particle.Trail, trailStart, trailEnd);
            }
        }

        // Draw particles
        foreach (var particle in allParticles)
        {
            Color4 color = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);

            // Add glow effect
            DrawCircle(particle.Body.Position, particle.Body.Shape is CircleShape c ? c.Radius * 2 : 0.3f,
                new Color4(color.R, color.G, color.B, 0.2f), true);

            // Draw particle
            DrawCircle(particle.Body.Position, particle.Body.Shape is CircleShape cs ? cs.Radius : 0.15f, color, true);

            // Draw velocity vector
            if (particle.Body.Velocity.Length > 0.5f)
            {
                ArtemisEngine.Vector2 velocityEnd = particle.Body.Position + particle.Body.Velocity.Normalized * 1.5f;
                DrawLine(particle.Body.Position, velocityEnd, new Color4(color.R, color.G, color.B, 0.5f), 1f);
            }
        }

        // Draw collision events
        foreach (var (evt, timestamp) in _recentCollisions)
        {
            float age = _experiment.SimulationTime - timestamp;
            if (age > 2.0f) continue;

            float alpha = 1f - (age / 2f);
            float radius = age * 4f + 0.5f;

            // Explosion effect
            Color4 explosionColor = new Color4(1f, 0.8f, 0.2f, alpha * 0.8f);
            DrawCircle(evt.CollisionPoint, radius, explosionColor, false);

            // Inner glow
            if (age < 0.5f)
            {
                DrawCircle(evt.CollisionPoint, radius * 0.5f, new Color4(1f, 1f, 1f, alpha), true);
            }
        }
    }

    protected override void DrawUI()
    {
        // UI is drawn in the title bar and console for now
        var (b1, b2, avgE) = _accelerator.GetStatistics();
        var particleCount = World.Bodies.Count(b => b.UserData is Particle);

        Title = $"Artemis Particle Collider | " +
                $"Beam1: {b1} | Beam2: {b2} | " +
                $"Energy: {avgE:F0}GeV | " +
                $"Collisions: {_experiment.TotalCollisions} | " +
                $"Particles: {particleCount} | " +
                $"TimeScale: {TimeScale:F2}x | " +
                (IsPaused ? "[PAUSED] " : "") +
                "[SPACE]Inject [1-3]Types [C]Clear [P]Pause [+-]Speed [ESC]Quit";
    }
}
