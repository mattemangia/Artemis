using ArtemisEngine;

namespace ParticleColliderDemo;

class Program
{
    static void Main(string[] args)
    {
        Console.Clear();
        Console.CursorVisible = false;
        Console.WriteLine("╔══════════════════════════════════════════════════════════════════════════════╗");
        Console.WriteLine("║                    ARTEMIS PARTICLE COLLIDER SIMULATOR                       ║");
        Console.WriteLine("║                        LHC-Style Event Display                               ║");
        Console.WriteLine("╚══════════════════════════════════════════════════════════════════════════════╝");
        Console.WriteLine();
        Console.WriteLine("This simulation demonstrates the scientific physics capabilities of Artemis:");
        Console.WriteLine("  • Continuous Collision Detection (CCD) for high-speed particles");
        Console.WriteLine("  • Energy tracking and conservation monitoring");
        Console.WriteLine("  • Advanced integration methods (RK4)");
        Console.WriteLine("  • Data export for scientific analysis");
        Console.WriteLine("  • Deterministic physics for reproducible experiments");
        Console.WriteLine();
        Console.WriteLine("Press any key to start the simulation...");
        Console.ReadKey(true);

        RunSimulation();
    }

    static void RunSimulation()
    {
        // Setup
        const int WIDTH = 80;
        const int HEIGHT = 40;
        const float DELTA_TIME = 1 / 60f;

        // Create physics world (no gravity - we're in vacuum!)
        var world = new PhysicsWorld(new Vector2(0, 0));
        world.UseAdvancedSolver = true; // Use Sequential Impulse solver

        // Create accelerator (27 units radius, like LHC's 27km)
        var accelerator = new Accelerator(
            center: new Vector2(0, 0),
            radius: 20f,
            beamEnergy: 100f // 100 GeV
        );

        // Create experiment
        var experiment = new CollisionExperiment("ARTEMIS Particle Physics Experiment", accelerator);

        // Create renderer
        var renderer = new LHCRenderer(WIDTH, HEIGHT);
        renderer.SetCamera(new Vector2(0, 0), zoom: 1.5f);

        // Add magnetic field to world
        world.AddAreaEffector(accelerator.GetMagneticField());

        // Add detectors to world
        foreach (var detector in experiment.Detectors)
        {
            world.AddBody(detector.TriggerZone);
        }

        // Setup collision listener for detector triggers
        world.OnTriggerEnter += (sender, e) =>
        {
            if (e.BodyA.UserData is Detector detector && e.BodyB.UserData is Particle particle)
            {
                detector.RecordDetection(particle, experiment.SimulationTime);
            }
            else if (e.BodyB.UserData is Detector detector2 && e.BodyA.UserData is Particle particle2)
            {
                detector2.RecordDetection(particle2, experiment.SimulationTime);
            }
        };

        // Track collision events
        List<(ParticleCollisionEvent evt, float timestamp)> recentCollisions = new();

        world.OnCollisionEnter += (sender, e) =>
        {
            if (e.BodyA.UserData is Particle p1 && e.BodyB.UserData is Particle p2)
            {
                Vector2 collisionPoint = (p1.Body.Position + p2.Body.Position) * 0.5f;
                experiment.RecordCollision(p1, p2, collisionPoint);

                var evt = new ParticleCollisionEvent(p1, p2, collisionPoint, experiment.SimulationTime);
                recentCollisions.Add((evt, experiment.SimulationTime));

                // Create collision products (simplified)
                if (Random.Shared.NextSingle() < 0.3f) // 30% chance
                {
                    CreateCollisionProducts(world, collisionPoint, p1, p2);
                }
            }
        };

        // Initial beam injection
        InjectInitialBeams(world, accelerator);

        // Main simulation loop
        bool running = true;
        DateTime lastFrame = DateTime.Now;

        while (running)
        {
            // Input handling
            if (Console.KeyAvailable)
            {
                var key = Console.ReadKey(true).Key;
                switch (key)
                {
                    case ConsoleKey.Q:
                        running = false;
                        break;

                    case ConsoleKey.Spacebar:
                        InjectBeams(world, accelerator);
                        break;

                    case ConsoleKey.C:
                        ClearExperiment(world, accelerator, experiment);
                        break;

                    case ConsoleKey.E:
                        ExportData(experiment);
                        break;

                    case ConsoleKey.D1:
                        InjectSpecificBeam(world, accelerator, ParticleType.Proton, ParticleType.Proton);
                        break;

                    case ConsoleKey.D2:
                        InjectSpecificBeam(world, accelerator, ParticleType.Electron, ParticleType.Positron);
                        break;

                    case ConsoleKey.D3:
                        InjectSpecificBeam(world, accelerator, ParticleType.Proton, ParticleType.Electron);
                        break;
                }
            }

            // Update physics
            world.Step(DELTA_TIME);

            // Update particles
            var allParticles = world.Bodies
                .Where(b => b.UserData is Particle)
                .Select(b => (Particle)b.UserData!)
                .ToList();

            foreach (var particle in allParticles)
            {
                particle.Update(DELTA_TIME);

                // Apply beam focusing
                accelerator.ApplyBeamFocusing(particle, DELTA_TIME);

                // Remove decayed particles
                if (particle.HasDecayed)
                {
                    world.RemoveBody(particle.Body);

                    // Add decay products
                    foreach (var product in particle.DecayProducts)
                    {
                        world.AddBody(product.Body);
                    }
                }
            }

            // Update experiment
            experiment.Update(DELTA_TIME, world);

            // Clean up old collision events (keep only last 2 seconds)
            recentCollisions.RemoveAll(c => experiment.SimulationTime - c.timestamp > 2.0f);

            // Rendering
            renderer.Clear();
            renderer.DrawAcceleratorRing(accelerator);
            renderer.DrawDetectors(experiment.Detectors);

            // Draw particles
            foreach (var particle in allParticles)
            {
                renderer.DrawParticle(particle);
            }

            // Draw recent collision events
            foreach (var (evt, timestamp) in recentCollisions)
            {
                float age = experiment.SimulationTime - timestamp;
                renderer.DrawCollisionEvent(evt, age);
            }

            renderer.Render();
            renderer.DrawUI(experiment, allParticles.Count);

            // Frame rate control
            var elapsed = (DateTime.Now - lastFrame).TotalMilliseconds;
            if (elapsed < 16) // Cap at ~60 FPS
            {
                Thread.Sleep((int)(16 - elapsed));
            }
            lastFrame = DateTime.Now;
        }

        Console.Clear();
        Console.CursorVisible = true;
        Console.WriteLine("Simulation ended.");
        Console.WriteLine($"Total collisions: {experiment.TotalCollisions}");
        Console.WriteLine($"Total simulation time: {experiment.SimulationTime:F1} seconds");
    }

    static void InjectInitialBeams(PhysicsWorld world, Accelerator accelerator)
    {
        // Inject proton beams (like LHC)
        accelerator.InjectBeam1(ParticleType.Proton, 8);
        accelerator.InjectBeam2(ParticleType.Proton, 8);

        // Add all particles to physics world
        foreach (var particle in accelerator.GetAllParticles())
        {
            world.AddBody(particle.Body);
        }
    }

    static void InjectBeams(PhysicsWorld world, Accelerator accelerator)
    {
        // Random particle types
        ParticleType[] types = { ParticleType.Proton, ParticleType.Electron, ParticleType.Positron };
        var type1 = types[Random.Shared.Next(types.Length)];
        var type2 = types[Random.Shared.Next(types.Length)];

        accelerator.InjectBeam1(type1, 5);
        accelerator.InjectBeam2(type2, 5);

        foreach (var particle in accelerator.Beam1.TakeLast(5).Concat(accelerator.Beam2.TakeLast(5)))
        {
            world.AddBody(particle.Body);
        }
    }

    static void InjectSpecificBeam(PhysicsWorld world, Accelerator accelerator, ParticleType type1, ParticleType type2)
    {
        accelerator.InjectBeam1(type1, 6);
        accelerator.InjectBeam2(type2, 6);

        foreach (var particle in accelerator.Beam1.TakeLast(6).Concat(accelerator.Beam2.TakeLast(6)))
        {
            world.AddBody(particle.Body);
        }
    }

    static void CreateCollisionProducts(PhysicsWorld world, Vector2 point, Particle p1, Particle p2)
    {
        // Energy available for product creation
        float availableEnergy = p1.GetKineticEnergy() + p2.GetKineticEnergy();

        // Create 2-3 new particles
        int productCount = Random.Shared.Next(2, 4);

        for (int i = 0; i < productCount; i++)
        {
            // Random direction
            float angle = Random.Shared.NextSingle() * MathF.PI * 2;
            Vector2 direction = new Vector2(MathF.Cos(angle), MathF.Sin(angle));

            // Random energy distribution
            float energyFraction = Random.Shared.NextSingle() * 0.4f;
            float speed = MathF.Sqrt(2 * availableEnergy * energyFraction);

            // Create photon or other light particle
            ParticleType productType = Random.Shared.NextSingle() < 0.5f
                ? ParticleType.Photon
                : ParticleType.Electron;

            var product = new Particle(productType, point, direction * speed);
            world.AddBody(product.Body);
        }
    }

    static void ClearExperiment(PhysicsWorld world, Accelerator accelerator, CollisionExperiment experiment)
    {
        // Remove all particles
        var particleBodies = world.Bodies
            .Where(b => b.UserData is Particle)
            .ToList();

        foreach (var body in particleBodies)
        {
            world.RemoveBody(body);
        }

        accelerator.ClearBeams();

        // Clear detectors
        foreach (var detector in experiment.Detectors)
        {
            detector.Clear();
        }
    }

    static void ExportData(CollisionExperiment experiment)
    {
        try
        {
            string filename = $"particle_collision_data_{DateTime.Now:yyyyMMdd_HHmmss}.csv";
            experiment.ExportData(filename);

            Console.SetCursorPosition(0, 45);
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine($"Data exported to: {filename}");
            Console.ResetColor();
            Thread.Sleep(1000);
        }
        catch (Exception ex)
        {
            Console.SetCursorPosition(0, 45);
            Console.ForegroundColor = ConsoleColor.Red;
            Console.WriteLine($"Export failed: {ex.Message}");
            Console.ResetColor();
            Thread.Sleep(1000);
        }
    }
}
