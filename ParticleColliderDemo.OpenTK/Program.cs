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

    // Cavalier 3D perspective settings
    private const float CavalierAngle = 30f * MathF.PI / 180f; // 30 degree tilt
    private const float DepthScale = 0.5f; // Depth foreshortening
    private const float RingHeight = 3f; // Height of the accelerator tube
    private float _viewRotation = 0f; // Rotation around Y axis for view

    // Mouse control
    private Vector2 _lastMousePos;
    private bool _isRotating = false;
    private bool _isPanning = false;

    public ParticleColliderWindow(int width, int height, string title)
        : base(width, height, title)
    {
        // Slow down for visibility
        TimeScale = 0.3f;
    }

    /// <summary>
    /// Transform 2D physics coordinates to cavalier 3D projection
    /// </summary>
    private ArtemisEngine.Vector2 ToCavalier(ArtemisEngine.Vector2 pos, float height = 0)
    {
        // Rotate around view axis
        float cos = MathF.Cos(_viewRotation);
        float sin = MathF.Sin(_viewRotation);
        float rx = pos.X * cos - pos.Y * sin;
        float ry = pos.X * sin + pos.Y * cos;

        // Apply cavalier projection:
        // X stays as X
        // Y is compressed and offset by depth
        float projX = rx;
        float projY = ry * MathF.Cos(CavalierAngle) + height;

        return new ArtemisEngine.Vector2(projX, projY);
    }

    /// <summary>
    /// Get depth value for sorting (higher Y = further back)
    /// </summary>
    private float GetDepth(ArtemisEngine.Vector2 pos)
    {
        float cos = MathF.Cos(_viewRotation);
        float sin = MathF.Sin(_viewRotation);
        return pos.X * sin + pos.Y * cos;
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
        // === MOUSE CONTROLS ===
        var mousePos = MouseState.Position;

        // Right-click drag to rotate view
        if (MouseState.IsButtonDown(MouseButton.Right))
        {
            if (!_isRotating)
            {
                _isRotating = true;
                _lastMousePos = mousePos;
            }
            else
            {
                float deltaX = mousePos.X - _lastMousePos.X;
                _viewRotation -= deltaX * 0.005f;
                _lastMousePos = mousePos;
            }
        }
        else
        {
            _isRotating = false;
        }

        // Middle-click drag to pan
        if (MouseState.IsButtonDown(MouseButton.Middle))
        {
            if (!_isPanning)
            {
                _isPanning = true;
                _lastMousePos = mousePos;
            }
            else
            {
                float deltaX = mousePos.X - _lastMousePos.X;
                float deltaY = mousePos.Y - _lastMousePos.Y;
                CameraPosition = new Vector2d(
                    CameraPosition.X - deltaX * 0.1,
                    CameraPosition.Y + deltaY * 0.1
                );
                _lastMousePos = mousePos;
            }
        }
        else
        {
            _isPanning = false;
        }

        // Mouse wheel to zoom
        if (MouseState.ScrollDelta.Y != 0)
        {
            CameraZoom *= 1f - MouseState.ScrollDelta.Y * 0.1f;
            CameraZoom = Math.Clamp(CameraZoom, 5, 100);
        }

        // === KEYBOARD CONTROLS ===

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

        // View rotation controls (Q/E) - keyboard fallback
        if (KeyboardState.IsKeyDown(Keys.Q))
        {
            _viewRotation += 0.02f;
        }
        if (KeyboardState.IsKeyDown(Keys.E))
        {
            _viewRotation -= 0.02f;
        }

        // Export data
        if (KeyboardState.IsKeyPressed(Keys.R))
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
        // === CAVALIER 3D PERSPECTIVE RENDERING ===

        // Draw ground plane grid (for depth reference)
        DrawGroundGrid();

        // Draw accelerator ring in 3D (back half first for depth sorting)
        Draw3DAcceleratorRing(drawBackHalf: true);

        // Draw detectors in 3D
        foreach (var detector in _experiment.Detectors)
        {
            Color4 detectorColor = detector.Type switch
            {
                DetectorType.Tracker => new Color4(0.2f, 0.5f, 0.8f, 0.4f),
                DetectorType.Calorimeter => new Color4(0.5f, 0.8f, 0.2f, 0.4f),
                DetectorType.MuonChamber => new Color4(0.8f, 0.5f, 0.2f, 0.4f),
                DetectorType.InteractionPoint => new Color4(1f, 0.2f, 0.2f, 0.6f),
                _ => new Color4(0.5f, 0.5f, 0.5f, 0.4f)
            };

            // Draw detector as 3D cylinder projection
            Draw3DDetector(detector.Position, detector.Radius, detectorColor);
        }

        // Draw interaction points in 3D
        foreach (var ip in _accelerator.GetInteractionPoints())
        {
            var projected = ToCavalier(ip, RingHeight);
            DrawCircle(projected, 1.0f, new Color4(0.8f, 0.2f, 0.2f, 1f), true);
            DrawCircle(projected, 1.5f, new Color4(1f, 0.3f, 0.3f, 0.4f), false);
        }

        // Get all particles sorted by depth
        var allParticles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .OrderBy(p => GetDepth(p.Body.Position))
            .ToList();

        // Draw particle trails in 3D
        foreach (var particle in allParticles)
        {
            if (particle.Trail.Count > 1)
            {
                Color4 particleColor = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);
                Draw3DTrail(particle.Trail, particleColor);
            }
        }

        // Draw particles in 3D with shadows
        foreach (var particle in allParticles)
        {
            Color4 color = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);
            float radius = particle.Body.Shape is CircleShape cs ? cs.Radius : 0.15f;

            // Draw shadow on ground
            var shadowPos = ToCavalier(particle.Body.Position, 0);
            DrawCircle(shadowPos, radius * 0.8f, new Color4(0f, 0f, 0f, 0.3f), true);

            // Draw particle at ring height
            var particlePos = ToCavalier(particle.Body.Position, RingHeight);

            // Glow effect
            DrawCircle(particlePos, radius * 2.5f, new Color4(color.R, color.G, color.B, 0.15f), true);
            DrawCircle(particlePos, radius * 1.5f, new Color4(color.R, color.G, color.B, 0.3f), true);

            // Main particle
            DrawCircle(particlePos, radius, color, true);

            // Velocity vector
            if (particle.Body.Velocity.Length > 0.5f)
            {
                var velEnd = ToCavalier(
                    particle.Body.Position + particle.Body.Velocity.Normalized * 2f,
                    RingHeight);
                DrawLine(particlePos, velEnd, new Color4(color.R, color.G, color.B, 0.6f), 1.5f);
            }
        }

        // Draw front half of ring (on top of particles)
        Draw3DAcceleratorRing(drawBackHalf: false);

        // Draw collision events in 3D
        foreach (var (evt, timestamp) in _recentCollisions)
        {
            float age = _experiment.SimulationTime - timestamp;
            if (age > 2.0f) continue;

            float alpha = 1f - (age / 2f);
            float radius = age * 5f + 0.5f;

            var collisionPos = ToCavalier(evt.CollisionPoint, RingHeight);

            // Multiple expanding rings for 3D explosion effect
            for (int i = 0; i < 3; i++)
            {
                float ringRadius = radius * (1f + i * 0.3f);
                float ringAlpha = alpha * (1f - i * 0.3f);
                DrawCircle(collisionPos, ringRadius, new Color4(1f, 0.8f - i * 0.2f, 0.2f, ringAlpha * 0.6f), false);
            }

            // Central flash
            if (age < 0.3f)
            {
                DrawCircle(collisionPos, radius * 0.3f, new Color4(1f, 1f, 1f, alpha), true);
            }
        }

        // Draw CERN-style data panels
        DrawCERNInterface();
    }

    /// <summary>
    /// Draw ground reference grid
    /// </summary>
    private void DrawGroundGrid()
    {
        Color4 gridColor = new Color4(0.15f, 0.15f, 0.2f, 0.5f);
        float gridSize = 5f;
        int gridLines = 8;

        for (int i = -gridLines; i <= gridLines; i++)
        {
            // Horizontal lines
            var start = ToCavalier(new ArtemisEngine.Vector2(-gridLines * gridSize, i * gridSize), 0);
            var end = ToCavalier(new ArtemisEngine.Vector2(gridLines * gridSize, i * gridSize), 0);
            DrawLine(start, end, gridColor, 1f);

            // Vertical lines
            start = ToCavalier(new ArtemisEngine.Vector2(i * gridSize, -gridLines * gridSize), 0);
            end = ToCavalier(new ArtemisEngine.Vector2(i * gridSize, gridLines * gridSize), 0);
            DrawLine(start, end, gridColor, 1f);
        }
    }

    /// <summary>
    /// Draw the accelerator ring as a 3D tube
    /// </summary>
    private void Draw3DAcceleratorRing(bool drawBackHalf)
    {
        int segments = 64;
        float tubeRadius = 0.8f;

        for (int i = 0; i < segments; i++)
        {
            float angle1 = 2 * MathF.PI * i / segments;
            float angle2 = 2 * MathF.PI * (i + 1) / segments;

            // Check if this segment is in front or back
            float midAngle = (angle1 + angle2) / 2;
            float depth = MathF.Sin(midAngle + _viewRotation);
            bool isBack = depth < 0;

            if (isBack != drawBackHalf)
                continue;

            // Ring position
            var p1 = new ArtemisEngine.Vector2(
                MathF.Cos(angle1) * _accelerator.Radius,
                MathF.Sin(angle1) * _accelerator.Radius);
            var p2 = new ArtemisEngine.Vector2(
                MathF.Cos(angle2) * _accelerator.Radius,
                MathF.Sin(angle2) * _accelerator.Radius);

            // Project to cavalier view
            var proj1Top = ToCavalier(p1, RingHeight + tubeRadius);
            var proj1Bot = ToCavalier(p1, RingHeight - tubeRadius);
            var proj2Top = ToCavalier(p2, RingHeight + tubeRadius);
            var proj2Bot = ToCavalier(p2, RingHeight - tubeRadius);

            // Color based on depth (darker = further)
            float brightness = drawBackHalf ? 0.25f : 0.5f;
            Color4 tubeColor = new Color4(0.3f * brightness, 0.4f * brightness, 0.6f * brightness, 0.9f);

            // Draw tube segment edges
            DrawLine(proj1Top, proj2Top, tubeColor, 2f);
            DrawLine(proj1Bot, proj2Bot, tubeColor, 2f);

            // Draw vertical connectors occasionally
            if (i % 4 == 0)
            {
                DrawLine(proj1Top, proj1Bot, tubeColor, 1f);
            }
        }

        // Draw beam line (glowing center)
        for (int i = 0; i < segments; i++)
        {
            float angle1 = 2 * MathF.PI * i / segments;
            float angle2 = 2 * MathF.PI * (i + 1) / segments;

            float depth = MathF.Sin((angle1 + angle2) / 2 + _viewRotation);
            bool isBack = depth < 0;
            if (isBack != drawBackHalf)
                continue;

            var p1 = new ArtemisEngine.Vector2(
                MathF.Cos(angle1) * _accelerator.Radius,
                MathF.Sin(angle1) * _accelerator.Radius);
            var p2 = new ArtemisEngine.Vector2(
                MathF.Cos(angle2) * _accelerator.Radius,
                MathF.Sin(angle2) * _accelerator.Radius);

            var proj1 = ToCavalier(p1, RingHeight);
            var proj2 = ToCavalier(p2, RingHeight);

            float brightness = drawBackHalf ? 0.4f : 0.8f;
            Color4 beamColor = new Color4(0.2f * brightness, 0.6f * brightness, 1f * brightness, 0.6f);
            DrawLine(proj1, proj2, beamColor, 1.5f);
        }
    }

    /// <summary>
    /// Draw detector as 3D cylinder projection
    /// </summary>
    private void Draw3DDetector(ArtemisEngine.Vector2 pos, float radius, Color4 color)
    {
        var topPos = ToCavalier(pos, RingHeight + 2);
        var botPos = ToCavalier(pos, RingHeight - 2);

        // Draw ellipse at top and bottom
        DrawEllipse3D(pos, radius, RingHeight + 2, color);
        DrawEllipse3D(pos, radius, RingHeight - 2, new Color4(color.R * 0.5f, color.G * 0.5f, color.B * 0.5f, color.A));

        // Draw vertical lines
        int vLines = 8;
        for (int i = 0; i < vLines; i++)
        {
            float angle = 2 * MathF.PI * i / vLines;
            var offset = new ArtemisEngine.Vector2(MathF.Cos(angle) * radius, MathF.Sin(angle) * radius);
            var top = ToCavalier(pos + offset, RingHeight + 2);
            var bot = ToCavalier(pos + offset, RingHeight - 2);
            DrawLine(top, bot, new Color4(color.R * 0.7f, color.G * 0.7f, color.B * 0.7f, color.A * 0.5f), 1f);
        }
    }

    /// <summary>
    /// Draw ellipse in 3D space
    /// </summary>
    private void DrawEllipse3D(ArtemisEngine.Vector2 center, float radius, float height, Color4 color)
    {
        int segments = 24;
        for (int i = 0; i < segments; i++)
        {
            float angle1 = 2 * MathF.PI * i / segments;
            float angle2 = 2 * MathF.PI * (i + 1) / segments;

            var p1 = center + new ArtemisEngine.Vector2(MathF.Cos(angle1) * radius, MathF.Sin(angle1) * radius);
            var p2 = center + new ArtemisEngine.Vector2(MathF.Cos(angle2) * radius, MathF.Sin(angle2) * radius);

            DrawLine(ToCavalier(p1, height), ToCavalier(p2, height), color, 1f);
        }
    }

    /// <summary>
    /// Draw particle trail in 3D
    /// </summary>
    private void Draw3DTrail(List<ArtemisEngine.Vector2> points, Color4 color)
    {
        if (points.Count < 2) return;

        for (int i = 0; i < points.Count - 1; i++)
        {
            float t = (float)i / points.Count;
            float alpha = t * 0.6f;
            Color4 trailColor = new Color4(color.R * t, color.G * t, color.B * t, alpha);

            var p1 = ToCavalier(points[i], RingHeight);
            var p2 = ToCavalier(points[i + 1], RingHeight);
            DrawLine(p1, p2, trailColor, 1f);
        }
    }

    /// <summary>
    /// Draw CERN-style interface panels (simulated with drawing primitives)
    /// </summary>
    private void DrawCERNInterface()
    {
        // We'll use the title bar for now, but this is where you'd add
        // overlay panels, graphs, and data displays in a full implementation
    }

    protected override void DrawUI()
    {
        // UI is drawn in the title bar
        var (b1, b2, avgE) = _accelerator.GetStatistics();
        var particleCount = World.Bodies.Count(b => b.UserData is Particle);

        Title = $"ARTEMIS LHC Simulator | " +
                $"Beam1: {b1} | Beam2: {b2} | " +
                $"E: {avgE:F0}GeV | " +
                $"Collisions: {_experiment.TotalCollisions} | " +
                $"Particles: {particleCount} | " +
                (IsPaused ? "[PAUSED] " : "") +
                "[SPACE]Inject [RightDrag]Rotate [Scroll]Zoom [C]Clear";
    }
}
