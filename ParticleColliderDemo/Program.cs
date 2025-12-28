using Artemis.Physics2D;
using Artemis.Physics2D.Particles;

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
        Console.WriteLine("  • Multi-threaded parallel physics simulation");
        Console.WriteLine("  • Continuous Collision Detection (CCD) for high-speed particles");
        Console.WriteLine("  • Energy tracking and conservation monitoring");
        Console.WriteLine("  • SIMD-optimized integration");
        Console.WriteLine("  • Data export for scientific analysis");
        Console.WriteLine();
        Console.WriteLine("Starting simulation in 3 seconds...");
        Thread.Sleep(3000);

        RunSimulation();
    }

    static void RunSimulation()
    {
        // Setup - use smaller viewport to prevent console scrolling
        const int WIDTH = 80;
        const int HEIGHT = 24;
        const double DELTA_TIME = 1.0 / 60.0;
        const double TIME_SCALE = 0.1; // Slow down simulation for visibility

        // Create physics world (no gravity - we're in vacuum!)
        var world = new PhysicsWorld2D(Vector2D.Zero);
        world.UseCCD = true; // Enable CCD for high-speed particles

        // Create accelerator (27 units radius, like LHC's 27km)
        var accelerator = new Accelerator(
            center: Vector2D.Zero,
            radius: 12.0,
            beamEnergy: 5.0
        );

        // Create experiment
        var experiment = new CollisionExperiment("ARTEMIS Particle Physics Experiment", accelerator);

        // Create renderer
        var renderer = new LHCRenderer(WIDTH, HEIGHT);
        renderer.SetCamera(Vector2D.Zero, zoom: 0.8);

        // Add magnetic field to world - Handled in Accelerator.Update via Lorentz force
        // world.AddAreaEffector(accelerator.GetMagneticField());

        // Add detectors to world
        foreach (var detector in experiment.Detectors)
        {
            if (detector.TriggerZone != null)
                world.AddBody(detector.TriggerZone);
        }

        // Setup collision listener for detector triggers
        world.TriggerEnter += (sender, e) =>
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
        List<(ParticleCollisionEvent evt, double timestamp)> recentCollisions = new();
        System.Collections.Concurrent.ConcurrentQueue<Particle> pendingParticles = new();

        world.CollisionBegin += (sender, e) =>
        {
            if (e.BodyA.UserData is Particle p1 && e.BodyB.UserData is Particle p2)
            {
                Vector2D collisionPoint = (p1.Body.Position + p2.Body.Position) * 0.5;
                experiment.RecordCollision(p1, p2, collisionPoint);

                var evt = new ParticleCollisionEvent(p1, p2, collisionPoint, experiment.SimulationTime);
                lock (recentCollisions)
                {
                    recentCollisions.Add((evt, experiment.SimulationTime));
                }

                // Create collision products (simplified)
                if (Random.Shared.NextDouble() < 0.3) // 30% chance
                {
                    var newParticles = CreateCollisionProducts(collisionPoint, p1, p2);
                    foreach (var p in newParticles)
                    {
                        pendingParticles.Enqueue(p);
                    }
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

            // Update physics with time scaling for visibility
            world.Step(DELTA_TIME * TIME_SCALE);

            while (pendingParticles.TryDequeue(out var p))
            {
                world.AddBody(p.Body);
            }

            // Update particles
            var allParticles = world.Bodies
                .Where(b => b.UserData is Particle)
                .Select(b => (Particle)b.UserData!)
                .ToList();

            // Update accelerator physics (acceleration + B-field)
            accelerator.Update(DELTA_TIME);

            foreach (var particle in allParticles)
            {
                particle.Update(DELTA_TIME);

                // Apply beam focusing
                accelerator.ApplyBeamFocusing(particle, DELTA_TIME);

                // Remove decayed particles
                if (particle.HasDecayed)
                {
                    world.RemoveBody(particle.Body);

                    // Add decay products (they use Particle2D internally)
                    foreach (var product in particle.DecayProducts)
                    {
                        world.AddBody(product.Body);
                    }
                }
            }

            // Update experiment
            experiment.Update(DELTA_TIME, world);

            // Clean up old collision events (keep only last 2 seconds)
            recentCollisions.RemoveAll(c => experiment.SimulationTime - c.timestamp > 2.0);

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
                double age = experiment.SimulationTime - timestamp;
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

    static void InjectInitialBeams(PhysicsWorld2D world, Accelerator accelerator)
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

    static void InjectBeams(PhysicsWorld2D world, Accelerator accelerator)
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

    static void InjectSpecificBeam(PhysicsWorld2D world, Accelerator accelerator, ParticleType type1, ParticleType type2)
    {
        accelerator.InjectBeam1(type1, 6);
        accelerator.InjectBeam2(type2, 6);

        foreach (var particle in accelerator.Beam1.TakeLast(6).Concat(accelerator.Beam2.TakeLast(6)))
        {
            world.AddBody(particle.Body);
        }
    }

    static List<Particle> CreateCollisionProducts(Vector2D point, Particle p1, Particle p2)
    {
        var products = new List<Particle>();

        // Energy available for product creation
        double availableEnergy = p1.GetKineticEnergy() + p2.GetKineticEnergy();

        // Create 2-3 new particles
        int productCount = Random.Shared.Next(2, 4);

        for (int i = 0; i < productCount; i++)
        {
            // Random direction
            double angle = Random.Shared.NextDouble() * Math.PI * 2;
            Vector2D direction = Vector2D.FromAngle(angle);

            // Random energy distribution
            double energyFraction = Random.Shared.NextDouble() * 0.4;
            double speed = Math.Sqrt(2 * availableEnergy * energyFraction);

            // Create photon or other light particle
            ParticleType productType = Random.Shared.NextDouble() < 0.5
                ? ParticleType.Photon
                : ParticleType.Electron;

            var product = new Particle(productType, point, direction * speed);
            products.Add(product);
        }

        return products;
    }

    static void ClearExperiment(PhysicsWorld2D world, Accelerator accelerator, CollisionExperiment experiment)
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
