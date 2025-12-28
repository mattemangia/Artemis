using Artemis.Graphics;
using Artemis.Physics2D;
using Artemis.Physics2D.Particles;
using OpenTK.Mathematics;
using OpenTK.Windowing.GraphicsLibraryFramework;
using ParticleColliderDemo;
using Vector2 = OpenTK.Mathematics.Vector2;

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
    private Vector2 ToCavalier(Vector2D pos, float height = 0)
    {
        // Rotate around view axis
        float cos = MathF.Cos(_viewRotation);
        float sin = MathF.Sin(_viewRotation);
        double rx = pos.X * cos - pos.Y * sin;
        double ry = pos.X * sin + pos.Y * cos;

        // Apply cavalier projection:
        // X stays as X
        // Y is compressed and offset by depth
        double projX = rx;
        double projY = ry * MathF.Cos(CavalierAngle) + height;

        return new Vector2(projX, projY);
    }

    /// <summary>
    /// Get depth value for sorting (higher Y = further back)
    /// </summary>
    private double GetDepth(Vector2D pos)
    {
        float cos = MathF.Cos(_viewRotation);
        float sin = MathF.Sin(_viewRotation);
        return pos.X * sin + pos.Y * cos;
    }

    protected override void Initialize()
    {
        // No gravity in the accelerator (vacuum)
        World = new PhysicsWorld2D(Vector2D.Zero);
        World.UseCCD = true;

        // Create accelerator with smaller energy for slower particles
        _accelerator = new Accelerator(
            center: Vector2D.Zero,
            radius: 15.0,
            beamEnergy: 3.0 // Lower energy = slower particles
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
        World.TriggerEnter += (sender, e) =>
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

        World.CollisionBegin += (sender, e) =>
        {
            if (e.BodyA.UserData is Particle p1 && e.BodyB.UserData is Particle p2)
            {
                Vector2D collisionPoint = (p1.Body.Position + p2.Body.Position) * 0.5;
                _experiment.RecordCollision(p1, p2, collisionPoint);

                var evt = new ParticleCollisionEvent(p1, p2, collisionPoint, _experiment.SimulationTime);
                _recentCollisions.Add((evt, (float)_experiment.SimulationTime));

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

    private void CreateCollisionProducts(Vector2D point, Particle p1, Particle p2)
    {
        // Limit total particles to prevent explosion
        int currentParticleCount = World.Bodies.Count(b => b.UserData is Particle);
        if (currentParticleCount > 100)
            return; // Don't create more if we already have too many

        double availableEnergy = p1.GetKineticEnergy() + p2.GetKineticEnergy();
        int productCount = Random.Shared.Next(2, 3); // Reduced from 2-4 to 2-3

        for (int i = 0; i < productCount; i++)
        {
            double angle = Random.Shared.NextDouble() * Math.PI * 2;
            Vector2D direction = Vector2D.FromAngle(angle);

            double energyFraction = Random.Shared.NextDouble() * 0.4;
            double speed = Math.Sqrt(2 * availableEnergy * energyFraction);

            ParticleType productType = Random.Shared.NextDouble() < 0.5
                ? ParticleType.Photon
                : ParticleType.Electron;

            // Spawn particles slightly away from collision point to prevent immediate re-collision
            Vector2D spawnPoint = point + direction * 0.5;
            var product = new Particle(productType, spawnPoint, direction * speed);

            // Disable CCD for collision products - they're slower and don't need it
            // This also prevents the infinite collision loop
            product.Body.IsBullet = false;

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
        _recentCollisions.RemoveAll(c => _experiment.SimulationTime - c.timestamp > 2.0);
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
                DetectorType2D.Tracker => new Color4(0.2f, 0.5f, 0.8f, 0.4f),
                DetectorType2D.Calorimeter => new Color4(0.5f, 0.8f, 0.2f, 0.4f),
                DetectorType2D.MuonChamber => new Color4(0.8f, 0.5f, 0.2f, 0.4f),
                DetectorType2D.InteractionPoint => new Color4(1f, 0.2f, 0.2f, 0.6f),
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

        // Draw particles as CERN-style point markers with velocity lines
        foreach (var particle in allParticles)
        {
            Color4 color = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);

            // Draw particle position as small point marker
            var particlePos = ToCavalier(particle.Body.Position, RingHeight);

            // Small point marker (like CERN event display)
            DrawCircle(particlePos, 0.12f, color, true);

            // Draw velocity/momentum vector as line
            double speed = particle.Body.Velocity.Magnitude;
            if (speed > 0.1)
            {
                // Line showing direction of movement
                var velEnd = ToCavalier(
                    particle.Body.Position + particle.Body.Velocity.Normalized * Math.Min(3.0, speed * 0.3),
                    RingHeight);
                DrawLine(particlePos, velEnd, color, 1.5f);
            }
        }

        // Draw front half of ring (on top of particles)
        Draw3DAcceleratorRing(drawBackHalf: false);

        // Draw collision events in 3D
        foreach (var (evt, timestamp) in _recentCollisions)
        {
            float age = (float)(_experiment.SimulationTime - timestamp);
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
            var start = ToCavalier(new Vector2(-gridLines * gridSize, i * gridSize), 0);
            var end = ToCavalier(new Vector2(gridLines * gridSize, i * gridSize), 0);
            DrawLine(start, end, gridColor, 1f);

            // Vertical lines
            start = ToCavalier(new Vector2(i * gridSize, -gridLines * gridSize), 0);
            end = ToCavalier(new Vector2(i * gridSize, gridLines * gridSize), 0);
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
            var p1 = new Vector2(
                MathF.Cos(angle1) * _accelerator.Radius,
                MathF.Sin(angle1) * _accelerator.Radius);
            var p2 = new Vector2(
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

            var p1 = new Vector2(
                MathF.Cos(angle1) * _accelerator.Radius,
                MathF.Sin(angle1) * _accelerator.Radius);
            var p2 = new Vector2(
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
    private void Draw3DDetector(Vector2D pos, double radius, Color4 color)
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
            double angle = 2 * Math.PI * i / vLines;
            var offset = new Vector2(Math.Cos(angle) * radius, Math.Sin(angle) * radius);
            var top = ToCavalier(pos + offset, RingHeight + 2);
            var bot = ToCavalier(pos + offset, RingHeight - 2);
            DrawLine(top, bot, new Color4(color.R * 0.7f, color.G * 0.7f, color.B * 0.7f, color.A * 0.5f), 1f);
        }
    }

    /// <summary>
    /// Draw ellipse in 3D space
    /// </summary>
    private void DrawEllipse3D(Vector2D center, double radius, float height, Color4 color)
    {
        int segments = 24;
        var centerV = ToVec2(center);
        for (int i = 0; i < segments; i++)
        {
            double angle1 = 2 * Math.PI * i / segments;
            double angle2 = 2 * Math.PI * (i + 1) / segments;

            var p1 = new Vector2D(centerV.X + Math.Cos(angle1) * radius, centerV.Y + Math.Sin(angle1) * radius);
            var p2 = new Vector2D(centerV.X + Math.Cos(angle2) * radius, centerV.Y + Math.Sin(angle2) * radius);

            DrawLine(ToCavalier(p1, height), ToCavalier(p2, height), color, 1f);
        }
    }

    /// <summary>
    /// Draw particle trail in 3D - CERN-style curved tracks
    /// </summary>
    private void Draw3DTrail(IList<Vector2D> points, Color4 color)
    {
        if (points.Count < 2) return;

        // Draw main track with fade
        for (int i = 0; i < points.Count - 1; i++)
        {
            float t = (float)i / points.Count;
            float alpha = 0.1f + t * 0.7f; // Start faint, get brighter
            float brightness = 0.3f + t * 0.7f;
            Color4 trailColor = new Color4(
                color.R * brightness,
                color.G * brightness,
                color.B * brightness,
                alpha);

            var p1 = ToCavalier(points[i], RingHeight);
            var p2 = ToCavalier(points[i + 1], RingHeight);
            DrawLine(p1, p2, trailColor, 1.5f);
        }

        // Add glow effect on recent part of trail
        int glowStart = Math.Max(0, points.Count - 10);
        for (int i = glowStart; i < points.Count - 1; i++)
        {
            float t = (float)(i - glowStart) / 10f;
            Color4 glowColor = new Color4(color.R, color.G, color.B, 0.1f + t * 0.2f);

            var p1 = ToCavalier(points[i], RingHeight);
            var p2 = ToCavalier(points[i + 1], RingHeight);
            DrawLine(p1, p2, glowColor, 3f);
        }
    }

    /// <summary>
    /// Draw CERN-style interface panels with text labels
    /// </summary>
    private void DrawCERNInterface()
    {
        // Calculate screen bounds in world coordinates for UI positioning
        float aspect = (float)Size.X / Size.Y;
        float halfWidth = (float)(CameraZoom * aspect);
        float halfHeight = (float)CameraZoom;

        float left = (float)CameraPosition.X - halfWidth;
        float right = (float)CameraPosition.X + halfWidth;
        float top = (float)CameraPosition.Y + halfHeight;
        float bottom = (float)CameraPosition.Y - halfHeight;

        // Panel dimensions
        float panelMargin = halfWidth * 0.02f;
        float panelWidth = halfWidth * 0.35f;
        float panelHeight = halfHeight * 0.5f;

        // Get statistics
        var (beam1Count, beam2Count, avgEnergy) = _accelerator.GetStatistics();
        int totalParticles = World.Bodies.Count(b => b.UserData is Particle);
        int totalCollisions = _experiment.TotalCollisions;

        // === LEFT PANEL - Beam Statistics ===
        float lpX = left + panelMargin;
        float lpY = top - panelMargin;

        // Panel background
        DrawPanelBackground(lpX, lpY, panelWidth, panelHeight);

        // Panel header
        float headerY = lpY - panelHeight * 0.08f;
        DrawLine(
            new Vector2(lpX + panelWidth * 0.05f, headerY),
            new Vector2(lpX + panelWidth * 0.95f, headerY),
            new Color4(0.2f, 0.8f, 1f, 0.8f), 2f);

        // BEAM 1 gauge with label
        float gaugeY = lpY - panelHeight * 0.2f;
        float gaugeWidth = panelWidth * 0.8f;
        float gaugeHeight = panelHeight * 0.08f;
        DrawGauge(lpX + panelWidth * 0.1f, gaugeY, gaugeWidth, gaugeHeight,
            beam1Count / 20f, new Color4(1f, 0.3f, 0.3f, 0.9f), "B1");

        // BEAM 2 gauge with label
        gaugeY -= panelHeight * 0.15f;
        DrawGauge(lpX + panelWidth * 0.1f, gaugeY, gaugeWidth, gaugeHeight,
            beam2Count / 20f, new Color4(0.3f, 0.6f, 1f, 0.9f), "B2");

        // Energy gauge with label
        gaugeY -= panelHeight * 0.2f;
        DrawGauge(lpX + panelWidth * 0.1f, gaugeY, gaugeWidth, gaugeHeight,
            MathF.Min(1f, (float)(avgEnergy / 1000f)), new Color4(1f, 0.8f, 0.2f, 0.9f), "E");

        // Particle count bars
        gaugeY -= panelHeight * 0.25f;
        DrawParticleCountBars(lpX + panelWidth * 0.1f, gaugeY, gaugeWidth, panelHeight * 0.2f);

        // === RIGHT PANEL - Collision Data ===
        float rpX = right - panelMargin - panelWidth;
        float rpY = top - panelMargin;

        // Panel background
        DrawPanelBackground(rpX, rpY, panelWidth, panelHeight * 0.7f);

        // Header line
        headerY = rpY - panelHeight * 0.08f;
        DrawLine(
            new Vector2(rpX + panelWidth * 0.05f, headerY),
            new Vector2(rpX + panelWidth * 0.95f, headerY),
            new Color4(1f, 0.5f, 0.2f, 0.8f), 2f);

        // Collision counter display (visual bars)
        float barY = rpY - panelHeight * 0.2f;
        int displayCollisions = Math.Min(totalCollisions, 50);
        DrawCollisionCounter(rpX + panelWidth * 0.1f, barY, panelWidth * 0.8f, panelHeight * 0.15f, displayCollisions);

        // Status indicators
        float statusY = rpY - panelHeight * 0.45f;
        DrawStatusIndicators(rpX + panelWidth * 0.1f, statusY, panelWidth * 0.8f);

        // === BOTTOM CENTER - ATLAS-style Event Display ===
        float ringCenterX = (float)CameraPosition.X;
        float ringCenterY = bottom + panelMargin + halfHeight * 0.2f;
        DrawATLASEventDisplay(ringCenterX, ringCenterY, halfWidth * 0.22f);

        // === TOP CENTER - Time Display ===
        DrawTimeDisplay((float)CameraPosition.X, top - panelMargin * 2, halfWidth * 0.3f);

        // === DRAW TEXT OVERLAYS (Screen Space) ===
        DrawTextOverlays(beam1Count, beam2Count, avgEnergy, totalParticles, totalCollisions);
    }

    /// <summary>
    /// Draw all text overlays for the CERN interface
    /// </summary>
    private void DrawTextOverlays(int beam1Count, int beam2Count, double avgEnergy, int totalParticles, int totalCollisions)
    {
        // Title
        DrawScreenText("ARTEMIS PARTICLE COLLIDER", 20, 20, 24, new Color4(0.3f, 0.8f, 1f, 1f));
        DrawScreenText("LHC Simulator", 20, 48, 14, new Color4(0.5f, 0.7f, 0.9f, 0.8f));

        // Left panel labels
        DrawScreenText("BEAM STATUS", 30, 85, 12, new Color4(0.2f, 0.8f, 1f, 0.9f));
        DrawScreenText($"BEAM 1: {beam1Count} particles", 40, 115, 11, new Color4(1f, 0.3f, 0.3f, 0.9f));
        DrawScreenText($"BEAM 2: {beam2Count} particles", 40, 155, 11, new Color4(0.3f, 0.6f, 1f, 0.9f));
        DrawScreenText($"ENERGY: {avgEnergy:F0} GeV", 40, 195, 11, new Color4(1f, 0.8f, 0.2f, 0.9f));
        DrawScreenText("PARTICLE COUNTS", 30, 235, 10, new Color4(0.6f, 0.8f, 1f, 0.7f));

        // Particle type legend
        string[] particleNames = { "p", "e-", "e+", "γ", "n", "q" };
        Color4[] particleColors = {
            new Color4(1f, 0.2f, 0.2f, 1f),  // Proton
            new Color4(0.2f, 0.4f, 1f, 1f),  // Electron
            new Color4(0.2f, 0.9f, 0.9f, 1f), // Positron
            new Color4(1f, 1f, 0.2f, 1f),    // Photon
            new Color4(0.6f, 0.6f, 0.6f, 1f), // Neutron
            new Color4(0.2f, 0.9f, 0.2f, 1f)  // Quark
        };
        for (int i = 0; i < particleNames.Length; i++)
        {
            DrawScreenText(particleNames[i], 45 + i * 30, 290, 10, particleColors[i]);
        }

        // Right panel labels
        int rightX = Size.X - 220;
        DrawScreenText("COLLISION DATA", rightX, 85, 12, new Color4(1f, 0.5f, 0.2f, 0.9f));
        DrawScreenText($"EVENTS: {totalCollisions}", rightX + 10, 110, 14, new Color4(1f, 0.8f, 0.3f, 1f));
        DrawScreenText($"PARTICLES: {totalParticles}", rightX + 10, 130, 11, new Color4(0.8f, 0.8f, 0.8f, 0.9f));

        // Status labels
        DrawScreenText("SYSTEM STATUS", rightX, 200, 10, new Color4(0.6f, 0.8f, 1f, 0.7f));
        string[] statusLabels = { "BEAM", "COLL", "HI-E", "DET", "SYS" };
        for (int i = 0; i < statusLabels.Length; i++)
        {
            DrawScreenText(statusLabels[i], rightX + 8 + i * 38, 245, 8, new Color4(0.5f, 0.6f, 0.7f, 0.8f));
        }

        // Bottom center - ATLAS display label
        int centerX = Size.X / 2;
        DrawScreenTextCentered("EVENT DISPLAY", centerX, Size.Y - 180, 11, new Color4(0.5f, 0.7f, 0.9f, 0.8f));
        DrawScreenTextCentered("ATLAS Cross-Section", centerX, Size.Y - 165, 9, new Color4(0.4f, 0.5f, 0.6f, 0.6f));

        // Time display
        float elapsed = (float)_experiment.SimulationTime;
        DrawScreenTextCentered($"T: {elapsed:F2}s", centerX, 60, 14, new Color4(0.2f, 0.6f, 1f, 0.9f));

        // Controls hint (bottom right)
        DrawScreenText("CONTROLS", Size.X - 200, Size.Y - 90, 10, new Color4(0.4f, 0.6f, 0.8f, 0.7f));
        DrawScreenText("[SPACE] Inject  [C] Clear", Size.X - 200, Size.Y - 75, 9, new Color4(0.5f, 0.5f, 0.5f, 0.7f));
        DrawScreenText("[1-3] Beam Types  [P] Pause", Size.X - 200, Size.Y - 62, 9, new Color4(0.5f, 0.5f, 0.5f, 0.7f));
        DrawScreenText("[RMB] Rotate  [Scroll] Zoom", Size.X - 200, Size.Y - 49, 9, new Color4(0.5f, 0.5f, 0.5f, 0.7f));
        DrawScreenText("[R] Export Data", Size.X - 200, Size.Y - 36, 9, new Color4(0.5f, 0.5f, 0.5f, 0.7f));

        // Collision event labels (in world space near collision points)
        foreach (var (evt, timestamp) in _recentCollisions)
        {
            float age = (float)(_experiment.SimulationTime - timestamp);
            if (age < 1.0f) // Show labels for recent collisions
            {
                var collisionPos = ToCavalier(evt.CollisionPoint, RingHeight);
                float alpha = 1f - age;
                var labelColor = new Color4(1f, 0.9f, 0.3f, alpha);

                double totalEnergy = evt.TotalEnergyBefore;
                DrawWorldText($"{totalEnergy:F0} GeV", new Vector2((float)collisionPos.X + 1f, (float)collisionPos.Y + 1f), 10, labelColor);
            }
        }
    }

    private void DrawPanelBackground(float x, float y, float width, float height)
    {
        // Dark semi-transparent background
        var center = new Vector2(x + width / 2, y - height / 2);
        DrawBox(center, width, height, 0, new Color4(0.05f, 0.08f, 0.12f, 0.85f), true);

        // Border
        DrawBox(center, width, height, 0, new Color4(0.2f, 0.4f, 0.6f, 0.6f), false);

        // Corner accents
        float cornerSize = width * 0.05f;
        Color4 accentColor = new Color4(0.3f, 0.7f, 1f, 0.8f);

        // Top-left corner
        DrawLine(new Vector2(x, y), new Vector2(x + cornerSize, y), accentColor, 2f);
        DrawLine(new Vector2(x, y), new Vector2(x, y - cornerSize), accentColor, 2f);

        // Top-right corner
        DrawLine(new Vector2(x + width, y), new Vector2(x + width - cornerSize, y), accentColor, 2f);
        DrawLine(new Vector2(x + width, y), new Vector2(x + width, y - cornerSize), accentColor, 2f);

        // Bottom-left corner
        DrawLine(new Vector2(x, y - height), new Vector2(x + cornerSize, y - height), accentColor, 2f);
        DrawLine(new Vector2(x, y - height), new Vector2(x, y - height + cornerSize), accentColor, 2f);

        // Bottom-right corner
        DrawLine(new Vector2(x + width, y - height), new Vector2(x + width - cornerSize, y - height), accentColor, 2f);
        DrawLine(new Vector2(x + width, y - height), new Vector2(x + width, y - height + cornerSize), accentColor, 2f);
    }

    private void DrawGauge(float x, float y, float width, float height, float fillPercent, Color4 fillColor, string label)
    {
        fillPercent = Math.Clamp(fillPercent, 0f, 1f);

        // Background
        var bgCenter = new Vector2(x + width / 2, y - height / 2);
        DrawBox(bgCenter, width, height, 0, new Color4(0.1f, 0.1f, 0.15f, 0.9f), true);

        // Fill bar
        if (fillPercent > 0.01f)
        {
            float fillWidth = width * fillPercent * 0.95f;
            var fillCenter = new Vector2(x + fillWidth / 2 + width * 0.025f, y - height / 2);
            DrawBox(fillCenter, fillWidth, height * 0.7f, 0, fillColor, true);

            // Glow effect on the fill
            DrawBox(fillCenter, fillWidth, height * 0.3f, 0, new Color4(1f, 1f, 1f, 0.3f), true);
        }

        // Border
        DrawBox(bgCenter, width, height, 0, new Color4(0.3f, 0.5f, 0.7f, 0.6f), false);

        // Label indicator (small box on left)
        var labelBox = new Vector2(x - height * 0.8f, y - height / 2);
        DrawBox(labelBox, height * 1.2f, height, 0, fillColor, true);
        DrawBox(labelBox, height * 1.2f, height, 0, new Color4(1f, 1f, 1f, 0.4f), false);
    }

    private void DrawParticleCountBars(float x, float y, float width, float height)
    {
        // Count particles by type
        var particles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .GroupBy(p => p.Properties.Type)
            .ToDictionary(g => g.Key, g => g.Count());

        ParticleType[] types = { ParticleType.Proton, ParticleType.Electron, ParticleType.Positron,
                                  ParticleType.Photon, ParticleType.Neutron, ParticleType.Quark };

        float barWidth = width / (types.Length + 1);
        float maxCount = 15f;

        for (int i = 0; i < types.Length; i++)
        {
            var type = types[i];
            int count = particles.GetValueOrDefault(type, 0);
            float barHeight = (count / maxCount) * height;
            barHeight = MathF.Min(barHeight, height);

            float barX = x + i * barWidth + barWidth * 0.5f;
            Color4 color = ParticleColors.GetValueOrDefault(type, Color4.White);

            if (count > 0)
            {
                var barCenter = new Vector2(barX, y - height + barHeight / 2);
                DrawBox(barCenter, barWidth * 0.7f, barHeight, 0, color, true);
            }

            // Base indicator
            DrawCircle(new Vector2(barX, y - height - barWidth * 0.3f), barWidth * 0.2f, color, true);
        }
    }

    private void DrawCollisionCounter(float x, float y, float width, float height, int count)
    {
        // Grid of collision indicator dots
        int cols = 10;
        int rows = 5;
        float dotRadius = MathF.Min(width / cols, height / rows) * 0.35f;

        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                int index = row * cols + col;
                float dotX = x + (col + 0.5f) * (width / cols);
                float dotY = y - (row + 0.5f) * (height / rows);

                Color4 dotColor;
                if (index < count)
                {
                    // Lit dot - gradient from yellow to red based on position
                    float t = (float)index / 50f;
                    dotColor = new Color4(1f, 1f - t * 0.6f, 0.2f - t * 0.2f, 0.9f);
                }
                else
                {
                    // Unlit dot
                    dotColor = new Color4(0.15f, 0.15f, 0.2f, 0.5f);
                }

                DrawCircle(new Vector2(dotX, dotY), dotRadius, dotColor, true);
            }
        }
    }

    private void DrawStatusIndicators(float x, float y, float width)
    {
        float indicatorSize = width * 0.08f;
        float spacing = width / 5f;

        // Status 1: Beam Active (green if particles exist)
        bool beamActive = World.Bodies.Any(b => b.UserData is Particle);
        DrawStatusLight(x + spacing * 0.5f, y, indicatorSize, beamActive,
            new Color4(0.2f, 1f, 0.3f, 1f), new Color4(0.1f, 0.3f, 0.1f, 0.5f));

        // Status 2: Collision Detected (orange if recent collision)
        bool recentCollision = _recentCollisions.Any(c => _experiment.SimulationTime - c.timestamp < 0.5f);
        DrawStatusLight(x + spacing * 1.5f, y, indicatorSize, recentCollision,
            new Color4(1f, 0.6f, 0.2f, 1f), new Color4(0.3f, 0.2f, 0.1f, 0.5f));

        // Status 3: High Energy (yellow if energy > 500)
        var (_, _, avgE) = _accelerator.GetStatistics();
        bool highEnergy = avgE > 500;
        DrawStatusLight(x + spacing * 2.5f, y, indicatorSize, highEnergy,
            new Color4(1f, 1f, 0.2f, 1f), new Color4(0.3f, 0.3f, 0.1f, 0.5f));

        // Status 4: Detector Active
        bool detectorActive = _experiment.Detectors.Any(d => d.DetectionCount > 0);
        DrawStatusLight(x + spacing * 3.5f, y, indicatorSize, detectorActive,
            new Color4(0.4f, 0.8f, 1f, 1f), new Color4(0.1f, 0.2f, 0.3f, 0.5f));

        // Status 5: System OK (always on)
        DrawStatusLight(x + spacing * 4.5f, y, indicatorSize, !IsPaused,
            new Color4(0.3f, 1f, 0.5f, 1f), new Color4(0.5f, 0.2f, 0.1f, 0.5f));
    }

    private void DrawStatusLight(float x, float y, float size, bool isOn, Color4 onColor, Color4 offColor)
    {
        Color4 color = isOn ? onColor : offColor;

        // Outer ring
        DrawCircle(new Vector2(x, y), size, new Color4(0.2f, 0.25f, 0.3f, 0.8f), true);

        // Inner light
        DrawCircle(new Vector2(x, y), size * 0.7f, color, true);

        // Highlight if on
        if (isOn)
        {
            DrawCircle(new Vector2(x, y), size * 1.3f, new Color4(color.R, color.G, color.B, 0.2f), true);
            DrawCircle(new Vector2(x - size * 0.2f, y + size * 0.2f), size * 0.2f,
                new Color4(1f, 1f, 1f, 0.4f), true);
        }
    }

    /// <summary>
    /// Draw CERN ATLAS-style event display with momentum/energy histogram
    /// </summary>
    private void DrawATLASEventDisplay(float centerX, float centerY, float radius)
    {
        var center = new Vector2(centerX, centerY);

        // Get particle statistics by angle sector
        var particles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .ToList();

        // Draw ATLAS-style calorimeter rings
        int sectors = 24;
        float[] sectorEnergies = new float[sectors];

        // Accumulate energy in each angular sector
        foreach (var particle in particles)
        {
            var relPos = particle.Body.Position - _accelerator.Center;
            double angle = Math.Atan2(relPos.Y, relPos.X);
            if (angle < 0) angle += 2 * Math.PI;
            int sector = (int)(angle / (2 * Math.PI) * sectors) % sectors;
            sectorEnergies[sector] += (float)particle.GetKineticEnergy();
        }

        // Find max energy for normalization
        float maxEnergy = sectorEnergies.Max();
        if (maxEnergy < 1) maxEnergy = 1;

        // Draw detector layers (like ATLAS cross-section)
        float innerRadius = radius * 0.3f;
        float midRadius = radius * 0.6f;
        float outerRadius = radius * 0.9f;

        // Inner tracker (blue)
        DrawRing(center, innerRadius, new Color4(0.1f, 0.3f, 0.6f, 0.6f), 2f);

        // Electromagnetic calorimeter (green)
        DrawRing(center, midRadius, new Color4(0.1f, 0.5f, 0.2f, 0.5f), 2f);

        // Hadronic calorimeter (yellow)
        DrawRing(center, outerRadius, new Color4(0.6f, 0.5f, 0.1f, 0.4f), 2f);

        // Draw energy deposits as sector bars (histogram style)
        for (int i = 0; i < sectors; i++)
        {
            float angle1 = (float)(2 * Math.PI * i / sectors) - MathF.PI / 2;
            float angle2 = (float)(2 * Math.PI * (i + 1) / sectors) - MathF.PI / 2;
            float midAngle = (angle1 + angle2) / 2;

            float energyRatio = sectorEnergies[i] / maxEnergy;
            if (energyRatio > 0.01f)
            {
                // Draw energy bar from inner to proportional outer radius
                float barLength = innerRadius + (outerRadius - innerRadius) * energyRatio;

                var innerPt = center + new Vector2(MathF.Cos(midAngle), MathF.Sin(midAngle)) * innerRadius;
                var outerPt = center + new Vector2(MathF.Cos(midAngle), MathF.Sin(midAngle)) * barLength;

                // Color based on energy (blue->green->yellow->red)
                Color4 barColor;
                if (energyRatio < 0.25f)
                    barColor = new Color4(0.2f, 0.4f, 1f, 0.8f);
                else if (energyRatio < 0.5f)
                    barColor = new Color4(0.2f, 0.8f, 0.4f, 0.8f);
                else if (energyRatio < 0.75f)
                    barColor = new Color4(0.9f, 0.8f, 0.2f, 0.8f);
                else
                    barColor = new Color4(1f, 0.3f, 0.2f, 0.9f);

                DrawLine(innerPt, outerPt, barColor, 3f);
            }
        }

        // Draw particle tracks through the detector (curved lines)
        foreach (var particle in particles.Take(20)) // Limit to prevent clutter
        {
            var relPos = particle.Body.Position - _accelerator.Center;
            double dist = relPos.Magnitude;

            // Only draw if particle is near the accelerator
            if (dist < _accelerator.Radius * 1.5)
            {
                double angle = Math.Atan2(relPos.Y, relPos.X);

                // Draw track from center to particle position
                var trackEnd = center + new Vector2(MathF.Cos((float)angle), MathF.Sin((float)angle)) *
                               (float)(innerRadius + (outerRadius - innerRadius) * Math.Min(1, dist / _accelerator.Radius));

                Color4 trackColor = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);
                trackColor = new Color4(trackColor.R, trackColor.G, trackColor.B, 0.4f);

                DrawLine(center, trackEnd, trackColor, 1f);
            }
        }

        // Center interaction point marker
        DrawCircle(center, radius * 0.08f, new Color4(1f, 0.2f, 0.2f, 0.9f), true);
        DrawCircle(center, radius * 0.12f, new Color4(1f, 0.4f, 0.4f, 0.4f), false);
    }

    private void DrawTimeDisplay(float centerX, float y, float width)
    {
        // Time elapsed bar
        float elapsed = (float)_experiment.SimulationTime;
        float cycleTime = elapsed % 10f; // 10 second cycles
        float cyclePercent = cycleTime / 10f;

        float barHeight = width * 0.08f;

        // Background
        var bgCenter = new Vector2(centerX, y - barHeight / 2);
        DrawBox(bgCenter, width, barHeight, 0, new Color4(0.08f, 0.1f, 0.15f, 0.8f), true);

        // Progress fill
        float fillWidth = width * cyclePercent * 0.95f;
        if (fillWidth > 0.1f)
        {
            var fillCenter = new Vector2(centerX - width / 2 + fillWidth / 2 + width * 0.025f, y - barHeight / 2);
            DrawBox(fillCenter, fillWidth, barHeight * 0.6f, 0, new Color4(0.2f, 0.6f, 1f, 0.7f), true);
        }

        // Tick marks
        for (int i = 0; i <= 10; i++)
        {
            float tickX = centerX - width / 2 + (width * i / 10f);
            float tickHeight = (i % 5 == 0) ? barHeight * 0.5f : barHeight * 0.25f;
            DrawLine(
                new Vector2(tickX, y),
                new Vector2(tickX, y + tickHeight),
                new Color4(0.4f, 0.6f, 0.8f, 0.6f), 1f);
        }

        // Border
        DrawBox(bgCenter, width, barHeight, 0, new Color4(0.2f, 0.4f, 0.6f, 0.5f), false);

        // Cycle indicator dots
        int cycleNumber = (int)(elapsed / 10f) % 5;
        for (int i = 0; i < 5; i++)
        {
            float dotX = centerX - width * 0.3f + i * width * 0.15f;
            float dotY = y + barHeight * 0.8f;
            bool lit = i <= cycleNumber;
            DrawCircle(new Vector2(dotX, dotY), barHeight * 0.2f,
                lit ? new Color4(0.3f, 0.8f, 1f, 0.9f) : new Color4(0.15f, 0.2f, 0.25f, 0.5f), true);
        }
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
