using Artemis.Graphics;
using Artemis.Physics2D;
using Artemis.Physics2D.Particles;
using OpenTK.Mathematics;
using OpenTK.Windowing.GraphicsLibraryFramework;
using OpenTK.Graphics.OpenGL4;
using ParticleColliderDemo;
// Using Artemis Vector2D by default for drawing as it seems GraphicsWindow expects it
using Vector2D = Artemis.Physics2D.Vector2D;

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

    // Beam history for graphs
    private readonly List<float> _beam1History = new();
    private readonly List<float> _beam2History = new();
    private const int MaxHistory = 100;
    private double _lastHistoryUpdate = 0;

    // Deferred actions to avoid LockRecursionException
    private readonly List<Action> _deferredActions = new();

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
    private const float RingHeight = 3f; // Height of the accelerator tube

    // View rotation
    private float _viewRotation = 0f; // Rotation around Y axis (Yaw)
    private float _viewPitch = 0f;    // Rotation around X axis (Pitch)

    // Mouse control
    private global::OpenTK.Mathematics.Vector2 _lastMousePos;
    private bool _isRotating = false;
    private bool _isPanning = false;

    public ParticleColliderWindow(int width, int height, string title)
        : base(width, height, title)
    {
        // Slow down for visibility
        TimeScale = 0.3f;
    }

    /// <summary>
    /// Transform 3D coordinates to 2D screen projection
    /// Uses rotation (Yaw) and Pitch
    /// </summary>
    private Vector2D ToPerspective(Vector2D pos, float height = 0)
    {
        // 1. Rotate around Y axis (Yaw)
        float cosY = MathF.Cos(_viewRotation);
        float sinY = MathF.Sin(_viewRotation);

        // 3D coordinates: X, Y (Height), Z (Depth)
        // Original 2D physics: X -> X, Y -> Z (Depth)
        float x = (float)pos.X;
        float z = (float)pos.Y; // Physics Y is Depth
        float y = height;       // Height is vertical

        // Rotate Y (Yaw)
        float rx = x * cosY - z * sinY;
        float rz = x * sinY + z * cosY;

        // 2. Rotate around X axis (Pitch)
        float cosP = MathF.Cos(_viewPitch);
        float sinP = MathF.Sin(_viewPitch);

        float rry = y * cosP - rz * sinP;
        float rrz = y * sinP + rz * cosP;

        // Simple perspective projection or orthographic with depth
        float scale = 1.0f;

        return new Vector2D(rx * scale, rry * scale);
    }

    // Overload for Vector2 input (OpenTK)
    private Vector2D ToPerspective(global::OpenTK.Mathematics.Vector2 pos, float height = 0)
        => ToPerspective(new Vector2D(pos.X, pos.Y), height);

    /// <summary>
    /// Get depth value for sorting (larger value = further away)
    /// </summary>
    private double GetDepth(Vector2D pos, float height = 0)
    {
        // Rotate around Y axis (Yaw)
        float cosY = MathF.Cos(_viewRotation);
        float sinY = MathF.Sin(_viewRotation);

        float x = (float)pos.X;
        float z = (float)pos.Y;
        float y = height;

        float rz = x * sinY + z * cosY;

        // Rotate around X axis (Pitch)
        float cosP = MathF.Cos(_viewPitch);
        float sinP = MathF.Sin(_viewPitch);

        float rrz = y * sinP + rz * cosP;

        return rrz;
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
            if (detector.TriggerZone != null)
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
                // Deferred to avoid lock recursion
                lock (_deferredActions)
                {
                    _deferredActions.Add(() =>
                    {
                        CreateCollisionProducts(collisionPoint, p1, p2);
                    });
                }
            }
        };

        // Inject initial beams
        InjectInitialBeams();

        // Set camera
        CameraZoom = 25;
        _viewPitch = 0.5f; // Initial slight tilt
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
        // Remove original particles (Consume them)
        World.RemoveBody(p1.Body);
        World.RemoveBody(p2.Body);

        _accelerator.Beam1.Remove(p1);
        _accelerator.Beam2.Remove(p1);
        _accelerator.Beam1.Remove(p2);
        _accelerator.Beam2.Remove(p2);

        // Limit total particles to prevent explosion
        int currentParticleCount = World.Bodies.Count(b => b.UserData is Particle);
        if (currentParticleCount > 100)
            return; // Don't create more if we already have too many

        double availableEnergy = p1.GetKineticEnergy() + p2.GetKineticEnergy();
        int productCount = Random.Shared.Next(2, 4); // 2-3 products

        for (int i = 0; i < productCount; i++)
        {
            double angle = Random.Shared.NextDouble() * Math.PI * 2;
            Vector2D direction = Vector2D.FromAngle(angle);

            double energyFraction = Random.Shared.NextDouble() * 0.4;
            double speed = Math.Sqrt(2 * availableEnergy * energyFraction);

            ParticleType productType = Random.Shared.NextDouble() < 0.5
                ? ParticleType.Photon
                : ParticleType.Electron;

            // Spawn particles slightly away from collision point
            Vector2D spawnPoint = point + direction * 0.5;
            var product = new Particle(productType, spawnPoint, direction * speed);

            // Disable CCD for collision products
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
                float deltaY = mousePos.Y - _lastMousePos.Y;

                _viewRotation -= deltaX * 0.005f; // Yaw
                _viewPitch -= deltaY * 0.005f;    // Pitch
                _viewPitch = Math.Clamp(_viewPitch, -MathF.PI / 2 + 0.1f, MathF.PI / 2 - 0.1f);

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
            lock (_deferredActions)
            {
                _deferredActions.Add(InjectBeams);
            }
        }

        // Clear experiment
        if (KeyboardState.IsKeyPressed(Keys.C))
        {
             lock (_deferredActions)
            {
                _deferredActions.Add(ClearExperiment);
            }
        }

        // Specific beam types
        if (KeyboardState.IsKeyPressed(Keys.D1))
        {
             lock (_deferredActions)
            {
                _deferredActions.Add(() => InjectSpecificBeam(ParticleType.Proton, ParticleType.Proton));
            }
        }
        if (KeyboardState.IsKeyPressed(Keys.D2))
        {
             lock (_deferredActions)
            {
                _deferredActions.Add(() => InjectSpecificBeam(ParticleType.Electron, ParticleType.Positron));
            }
        }
        if (KeyboardState.IsKeyPressed(Keys.D3))
        {
             lock (_deferredActions)
            {
                _deferredActions.Add(() => InjectSpecificBeam(ParticleType.Proton, ParticleType.Electron));
            }
        }

        // View rotation controls (Q/E/R/F)
        if (KeyboardState.IsKeyDown(Keys.Q)) _viewRotation += 0.02f;
        if (KeyboardState.IsKeyDown(Keys.E)) _viewRotation -= 0.02f;
        if (KeyboardState.IsKeyDown(Keys.R)) _viewPitch += 0.02f;
        if (KeyboardState.IsKeyDown(Keys.F)) _viewPitch -= 0.02f;

        // Export data
        if (KeyboardState.IsKeyPressed(Keys.T))
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
        // Execute deferred actions
        lock (_deferredActions)
        {
            foreach (var action in _deferredActions)
            {
                action();
            }
            _deferredActions.Clear();
        }

        // Update particles
        var allParticles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .ToList();

        // Update accelerator
        _accelerator.Update(deltaTime);

        foreach (var particle in allParticles)
        {
            particle.Update(deltaTime);

            // Apply beam focusing
            _accelerator.ApplyBeamFocusing(particle, deltaTime);

            // Remove decayed particles
            if (particle.HasDecayed)
            {
                lock (_deferredActions)
                {
                    _deferredActions.Add(() =>
                    {
                        if (World.RemoveBody(particle.Body))
                        {
                            foreach (var product in particle.DecayProducts)
                            {
                                World.AddBody(product.Body);
                            }
                        }
                    });
                }
            }
        }

        // Update experiment time
        _experiment.Update(deltaTime, World);

        // Update history for graphs (every 0.5s)
        if (_experiment.SimulationTime - _lastHistoryUpdate > 0.5)
        {
            var (b1, b2, _) = _accelerator.GetStatistics();
            _beam1History.Add(b1);
            _beam2History.Add(b2);

            if (_beam1History.Count > MaxHistory) _beam1History.RemoveAt(0);
            if (_beam2History.Count > MaxHistory) _beam2History.RemoveAt(0);

            _lastHistoryUpdate = _experiment.SimulationTime;
        }

        // Clean up old collision events
        _recentCollisions.RemoveAll(c => _experiment.SimulationTime - c.timestamp > 2.0);
    }

    protected override void Render()
    {
        // Split screen: Top 80% for 3D view, Bottom 20% for UI
        int uiHeight = (int)(Size.Y * 0.25f);
        int viewHeight = Size.Y - uiHeight;

        // === VIEWPORT 1: 3D Scene ===
        GL.Viewport(0, uiHeight, Size.X, viewHeight);

        float aspect = (float)Size.X / viewHeight;
        float halfWidth = (float)(CameraZoom * aspect);
        float halfHeight = (float)CameraZoom;

        // Re-calculate projection for top viewport
        Matrix4 projection = Matrix4.CreateOrthographicOffCenter(
            (float)CameraPosition.X - halfWidth,
            (float)CameraPosition.X + halfWidth,
            (float)CameraPosition.Y - halfHeight,
            (float)CameraPosition.Y + halfHeight,
            -1f, 1f);

        int projectionLoc = GL.GetUniformLocation(ShaderProgram, "projection");
        GL.UniformMatrix4(projectionLoc, false, ref projection);

        // Draw ground plane grid (for depth reference)
        DrawGroundGrid();

        // Draw accelerator ring in 3D (back half first for depth sorting)
        Draw3DAcceleratorRing(drawBackHalf: true);

        // Draw detectors (Back)
        foreach (var detector in _experiment.Detectors)
        {
            if (GetDepth(detector.Position, 0) < 0)
                Draw3DDetector(detector);
        }

        // Draw interaction points (Back)
        foreach (var ip in _accelerator.GetInteractionPoints())
        {
            if (GetDepth(ip, 0) < 0)
            {
                 var projected = ToPerspective(ip, RingHeight);
                 DrawCircle(projected, 1.0f, new Color4(0.8f, 0.2f, 0.2f, 1f), true);
            }
        }

        // Draw particles sorted by depth
        var allParticles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .OrderBy(p => GetDepth(p.Body.Position, RingHeight))
            .ToList();

        // Draw particle trails
        foreach (var particle in allParticles)
        {
            if (particle.Trail.Count > 1)
            {
                Color4 particleColor = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);
                Draw3DTrail(particle.Trail, particleColor);
            }
        }

        // Draw particles
        foreach (var particle in allParticles)
        {
            Color4 color = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);
            var particlePos = ToPerspective(particle.Body.Position, RingHeight);
            DrawCircle(particlePos, 0.12f, color, true);

            double speed = particle.Body.Velocity.Magnitude;
            if (speed > 0.1)
            {
                var endPos3D = particle.Body.Position + particle.Body.Velocity.Normalized * Math.Min(3.0, speed * 0.3);
                var velEnd = ToPerspective(endPos3D, RingHeight);
                DrawLine(particlePos, velEnd, color, 1.5f);
            }
        }

        // Draw front half of ring
        Draw3DAcceleratorRing(drawBackHalf: false);

        // Draw Front Detectors/IPs
        foreach (var detector in _experiment.Detectors)
        {
            if (GetDepth(detector.Position, 0) >= 0)
                Draw3DDetector(detector);
        }
         foreach (var ip in _accelerator.GetInteractionPoints())
        {
            if (GetDepth(ip, 0) >= 0)
            {
                 var projected = ToPerspective(ip, RingHeight);
                 DrawCircle(projected, 1.0f, new Color4(0.8f, 0.2f, 0.2f, 1f), true);
            }
        }

        // Draw collision events
        foreach (var (evt, timestamp) in _recentCollisions)
        {
            float age = (float)(_experiment.SimulationTime - timestamp);
            if (age > 2.0f) continue;

            float alpha = 1f - (age / 2f);
            float radius = age * 2f + 0.5f; // Reduced size

            var collisionPos = ToPerspective(evt.CollisionPoint, RingHeight);

            // Multiple expanding rings
            for (int i = 0; i < 3; i++)
            {
                float ringRadius = radius * (1f + i * 0.2f);
                float ringAlpha = alpha * (1f - i * 0.3f);
                DrawCircle(collisionPos, ringRadius, new Color4(1f, 0.8f - i * 0.2f, 0.2f, ringAlpha * 0.4f), false);
            }

            // Central flash
            if (age < 0.2f)
            {
                DrawCircle(collisionPos, radius * 0.4f, new Color4(1f, 1f, 1f, alpha * 0.8f), true);
            }
        }

        // Restore viewport to full window for any subsequent calls if any
        GL.Viewport(0, 0, Size.X, Size.Y);
    }

    /// <summary>
    /// Draw ground reference grid
    /// </summary>
    private void DrawGroundGrid()
    {
        Color4 gridColor = new Color4(0.15f, 0.15f, 0.2f, 0.5f);
        float gridSize = 5f;
        int gridLines = 8;

        // Draw grid lines
        for (int i = -gridLines; i <= gridLines; i++)
        {
            // Horizontal lines (along X)
            var start = ToPerspective(new Vector2D(-gridLines * gridSize, i * gridSize), -5);
            var end = ToPerspective(new Vector2D(gridLines * gridSize, i * gridSize), -5);
            DrawLine(start, end, gridColor, 1f);

            // Vertical lines (along Z)
            start = ToPerspective(new Vector2D(i * gridSize, -gridLines * gridSize), -5);
            end = ToPerspective(new Vector2D(i * gridSize, gridLines * gridSize), -5);
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

            float z1 = MathF.Sin(angle1) * (float)_accelerator.Radius;
            float z2 = MathF.Sin(angle2) * (float)_accelerator.Radius;
            float midX = MathF.Cos((angle1 + angle2)/2) * (float)_accelerator.Radius;
            float midZ = (z1 + z2) / 2;

            double depth = GetDepth(new Vector2D(midX, midZ), RingHeight);
            bool isBack = depth > 0;

            if (isBack != drawBackHalf)
                continue;

            // Ring position
            var p1 = new Vector2D(
                (float)(MathF.Cos(angle1) * _accelerator.Radius),
                (float)(MathF.Sin(angle1) * _accelerator.Radius));
            var p2 = new Vector2D(
                (float)(MathF.Cos(angle2) * _accelerator.Radius),
                (float)(MathF.Sin(angle2) * _accelerator.Radius));

            // Project to perspective view
            var proj1Top = ToPerspective(p1, RingHeight + tubeRadius);
            var proj1Bot = ToPerspective(p1, RingHeight - tubeRadius);
            var proj2Top = ToPerspective(p2, RingHeight + tubeRadius);
            var proj2Bot = ToPerspective(p2, RingHeight - tubeRadius);

            // Color based on depth (darker = further)
            float brightness = drawBackHalf ? 0.3f : 0.6f;
            Color4 tubeColor = new Color4(0.3f * brightness, 0.4f * brightness, 0.6f * brightness, 0.9f);

            // Draw tube segment edges
            DrawLine(proj1Top, proj2Top, tubeColor, 2f);
            DrawLine(proj1Bot, proj2Bot, tubeColor, 2f);

            // Draw sides
            DrawLine(proj1Top, proj1Bot, tubeColor, 1f);

            // Draw beam line (glowing center)
            var proj1 = ToPerspective(p1, RingHeight);
            var proj2 = ToPerspective(p2, RingHeight);

            float beamBrightness = drawBackHalf ? 0.4f : 0.8f;
            Color4 beamColor = new Color4(0.2f * beamBrightness, 0.6f * beamBrightness, 1f * beamBrightness, 0.6f);
            DrawLine(proj1, proj2, beamColor, 1.5f);
        }
    }

    /// <summary>
    /// Draw detector as 3D cylinder projection
    /// </summary>
    private void Draw3DDetector(Detector detector)
    {
        Vector2D pos = detector.Position;
        double radius = detector.Radius;

        Color4 color = detector.Type switch
            {
                DetectorType2D.Tracker => new Color4(0.2f, 0.5f, 0.8f, 0.4f),
                DetectorType2D.Calorimeter => new Color4(0.5f, 0.8f, 0.2f, 0.4f),
                DetectorType2D.MuonChamber => new Color4(0.8f, 0.5f, 0.2f, 0.4f),
                DetectorType2D.InteractionPoint => new Color4(1f, 0.2f, 0.2f, 0.6f),
                _ => new Color4(0.5f, 0.5f, 0.5f, 0.4f)
            };

        var topPos = ToPerspective(pos, RingHeight + 2);
        var botPos = ToPerspective(pos, RingHeight - 2);

        // Draw ellipse at top and bottom
        DrawEllipse3D(pos, radius, RingHeight + 2, color);
        DrawEllipse3D(pos, radius, RingHeight - 2, new Color4(color.R * 0.5f, color.G * 0.5f, color.B * 0.5f, color.A));

        // Draw vertical lines
        int vLines = 8;
        for (int i = 0; i < vLines; i++)
        {
            double angle = 2 * Math.PI * i / vLines;
            var offset = new Vector2D(Math.Cos(angle) * radius, Math.Sin(angle) * radius);
            // This offset is in XZ plane (Physics XY)
            var p = pos + offset;

            var top = ToPerspective(p, RingHeight + 2);
            var bot = ToPerspective(p, RingHeight - 2);
            DrawLine(top, bot, new Color4(color.R * 0.7f, color.G * 0.7f, color.B * 0.7f, color.A * 0.5f), 1f);
        }
    }

    /// <summary>
    /// Draw ellipse in 3D space (circle in XZ plane)
    /// </summary>
    private void DrawEllipse3D(Vector2D center, double radius, float height, Color4 color)
    {
        int segments = 24;
        for (int i = 0; i < segments; i++)
        {
            double angle1 = 2 * Math.PI * i / segments;
            double angle2 = 2 * Math.PI * (i + 1) / segments;

            var p1 = center + new Vector2D(Math.Cos(angle1) * radius, Math.Sin(angle1) * radius);
            var p2 = center + new Vector2D(Math.Cos(angle2) * radius, Math.Sin(angle2) * radius);

            DrawLine(ToPerspective(p1, height), ToPerspective(p2, height), color, 1f);
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

            var p1 = ToPerspective(points[i], RingHeight);
            var p2 = ToPerspective(points[i + 1], RingHeight);
            DrawLine(p1, p2, trailColor, 1.5f);
        }
    }

    protected override void DrawUI()
    {
        int uiHeight = (int)(Size.Y * 0.25f);
        GL.Viewport(0, 0, Size.X, uiHeight);

        // Orthographic projection for UI (0,0 bottom left to width,height top right)
        Matrix4 uiProjection = Matrix4.CreateOrthographicOffCenter(0, Size.X, 0, uiHeight, -1f, 1f);

        // Use a clean shader for UI
        int projectionLoc = GL.GetUniformLocation(ShaderProgram, "projection");
        GL.UniformMatrix4(projectionLoc, false, ref uiProjection);

        // Draw background for UI panel (Standard CERN Grey/Black)
        DrawBox(new Vector2D(Size.X/2, uiHeight/2), Size.X, uiHeight, 0, new Color4(0.85f, 0.85f, 0.85f, 1f), true); // Grey bg
        DrawBox(new Vector2D(Size.X/2, uiHeight/2), Size.X - 4, uiHeight - 4, 0, new Color4(0f, 0f, 0f, 1f), true);   // Black inner

        DrawCERNInterface(uiHeight);

        // Restore viewport
        GL.Viewport(0, 0, Size.X, Size.Y);

        // Set title
         var (b1, b2, avgE) = _accelerator.GetStatistics();
        Title = $"ARTEMIS LHC | B1:{b1} B2:{b2} | E:{avgE:F0}GeV | Collisions:{_experiment.TotalCollisions}";
    }

    /// <summary>
    /// Draw "LHC Page 1" style interface
    /// </summary>
    private void DrawCERNInterface(int height)
    {
        float width = Size.X;
        float margin = 10f;
        float w = width - 2 * margin;
        float h = height - 2 * margin;
        float startX = margin;
        float startY = margin;

        // Coordinates: Y is 0 at bottom in GL, increasing upwards.
        // ScreenText Y: 0 is Top.
        // We need to coordinate these.
        // Screen Text Y for bottom panel:
        float uiTopY = Size.Y - height + margin;

        // === HEADER ===
        // "LHC PAGE 1" - Top Left
        DrawScreenText("LHC PAGE 1", startX, uiTopY, 16, new Color4(0.2f, 0.2f, 1f, 1f));
        DrawScreenText($"Fill: {DateTime.Now.DayOfYear}{DateTime.Now.Hour}", startX + 150, uiTopY + 2, 14, Color4.White);

        var (b1, b2, avgE) = _accelerator.GetStatistics();
        DrawScreenText($"E: {avgE:F0} GeV", width - 150, uiTopY, 16, new Color4(0.2f, 1f, 1f, 1f));
        DrawScreenText(DateTime.Now.ToString("dd-MM-yy HH:mm:ss"), width - 220, uiTopY + 20, 12, new Color4(0.8f, 0.8f, 0.8f, 1f));

        // === STATUS ===
        // "STABLE BEAMS" - Center Top (Green background box with black text)
        string status = (b1 > 0 && b2 > 0) ? "STABLE BEAMS" : "INJECTION PHYSICS BEAM";
        Color4 statusColor = (b1 > 0 && b2 > 0) ? new Color4(0f, 0.8f, 0f, 1f) : new Color4(1f, 0.6f, 0f, 1f);

        float statusW = 300;
        float statusH = 30;
        float statusX = width / 2;
        float statusY = height - 30; // GL coordinates (near top of UI panel)

        DrawBox(new Vector2D(statusX, statusY), statusW, statusH, 0, statusColor, true);
        DrawScreenTextCentered(status, width / 2, uiTopY + 15, 14, Color4.Black);

        // === GRAPHS ===
        // Left: Beam 1 Intensity, Right: Beam 2 Intensity
        float graphW = w / 2 - 20;
        float graphH = h * 0.6f;
        float graphY = height / 2 - 10; // GL coords center

        float graph1X = width * 0.25f;
        float graph2X = width * 0.75f;

        // Draw Graph Frames (White border, Black bg)
        DrawBox(new Vector2D(graph1X, graphY), graphW, graphH, 0, new Color4(0.1f, 0.1f, 0.1f, 1f), true);
        DrawBox(new Vector2D(graph1X, graphY), graphW, graphH, 0, Color4.White, false);

        DrawBox(new Vector2D(graph2X, graphY), graphW, graphH, 0, new Color4(0.1f, 0.1f, 0.1f, 1f), true);
        DrawBox(new Vector2D(graph2X, graphY), graphW, graphH, 0, Color4.White, false);

        // Plot History
        DrawHistoryGraph(new Vector2D(graph1X, graphY), graphW, graphH, _beam1History, new Color4(0.2f, 0.2f, 1f, 1f));
        DrawHistoryGraph(new Vector2D(graph2X, graphY), graphW, graphH, _beam2History, new Color4(1f, 0.2f, 0.2f, 1f));

        // Labels under graphs
        // Screen text coordinates
        float labelY = uiTopY + height * 0.75f;
        DrawScreenText("I(B1):", graph1X - graphW/2 + 10, labelY, 12, Color4.White);
        DrawScreenText($"{b1:E2}", graph1X + 10, labelY, 12, new Color4(0.2f, 0.5f, 1f, 1f));

        DrawScreenText("I(B2):", graph2X - graphW/2 + 10, labelY, 12, Color4.White);
        DrawScreenText($"{b2:E2}", graph2X + 10, labelY, 12, new Color4(1f, 0.3f, 0.3f, 1f));

        // === COMMENTS ===
        DrawScreenText("Comments (02-Jan-2025 14:22:15):", startX, Size.Y - 30, 10, Color4.Gray);
        DrawScreenText("Collisions in IP1, IP2, IP5, IP8. Lumis optimization.", startX, Size.Y - 15, 12, new Color4(0.5f, 1f, 0.5f, 1f));

        // Event display overlay (small, center bottom)
        DrawATLASEventDisplay(width/2, graphY, height * 0.5f);
    }

    private void DrawHistoryGraph(Vector2D center, float w, float h, List<float> history, Color4 color)
    {
        if (history.Count < 2) return;

        float xStep = w / MaxHistory;
        float startX = (float)center.X - w/2;
        float bottomY = (float)center.Y - h/2;

        float maxVal = history.Max();
        if (maxVal < 1) maxVal = 1;

        for (int i = 0; i < history.Count - 1; i++)
        {
            float val1 = history[i];
            float val2 = history[i+1];

            float x1 = startX + i * xStep;
            float x2 = startX + (i+1) * xStep;

            float y1 = bottomY + (val1 / maxVal) * h * 0.9f;
            float y2 = bottomY + (val2 / maxVal) * h * 0.9f;

            DrawLine(new Vector2D(x1, y1), new Vector2D(x2, y2), color, 1.5f);
        }
    }

    private void DrawUIBar(float cx, float cy, float w, float h, float pct, Color4 color, string label)
    {
        pct = Math.Clamp(pct, 0, 1);
        DrawBox(new Vector2D(cx, cy), w, h, 0, new Color4(0.2f, 0.2f, 0.3f, 0.8f), true);
        if (pct > 0)
        {
            float fw = w * pct;
            DrawBox(new Vector2D(cx - w/2 + fw/2, cy), fw, h, 0, color, true);
        }
        DrawBox(new Vector2D(cx, cy), w, h, 0, Color4.Gray, false);
    }

    private void DrawCollisionCounter(float x, float y, float width, float height, int count)
    {
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

                Color4 dotColor = index < count
                    ? new Color4(1f, 1f - index/50f * 0.6f, 0.2f, 0.9f)
                    : new Color4(0.15f, 0.15f, 0.2f, 0.5f);

                DrawCircle(new Vector2D(dotX, dotY), dotRadius, dotColor, true);
            }
        }
    }

    private void DrawStatusIndicators(float x, float y, float width)
    {
        float indicatorSize = width * 0.08f;
        float spacing = width / 5f;

        bool beamActive = World.Bodies.Any(b => b.UserData is Particle);
        DrawStatusLight(x + spacing * 0.5f, y, indicatorSize, beamActive,
            new Color4(0.2f, 1f, 0.3f, 1f), new Color4(0.1f, 0.3f, 0.1f, 0.5f));

        bool recentCollision = _recentCollisions.Any(c => _experiment.SimulationTime - c.timestamp < 0.5f);
        DrawStatusLight(x + spacing * 1.5f, y, indicatorSize, recentCollision,
            new Color4(1f, 0.6f, 0.2f, 1f), new Color4(0.3f, 0.2f, 0.1f, 0.5f));
    }

    private void DrawStatusLight(float x, float y, float size, bool isOn, Color4 onColor, Color4 offColor)
    {
        Color4 color = isOn ? onColor : offColor;
        DrawCircle(new Vector2D(x, y), size, new Color4(0.2f, 0.25f, 0.3f, 0.8f), true);
        DrawCircle(new Vector2D(x, y), size * 0.7f, color, true);
    }

    private void DrawATLASEventDisplay(float centerX, float centerY, float radius)
    {
        var center = new Vector2D(centerX, centerY);

        int rings = 3;
        for (int i = 0; i < rings; i++)
        {
            float r = radius * (0.3f + i * 0.3f);
            float alpha = 0.5f - i * 0.1f;
            Color4 ringColor = i switch {
                0 => new Color4(0.2f, 0.4f, 1f, alpha),
                1 => new Color4(0.2f, 0.8f, 0.2f, alpha),
                _ => new Color4(0.8f, 0.6f, 0.2f, alpha)
            };
            DrawRing(center, r, ringColor, 2f);

            int spokes = 8;
            for (int s = 0; s < spokes; s++)
            {
                float angle = s * MathF.PI * 2 / spokes;
                var start = center + new Vector2D(MathF.Cos(angle), MathF.Sin(angle)) * r * 0.8f;
                var end = center + new Vector2D(MathF.Cos(angle), MathF.Sin(angle)) * r;
                DrawLine(start, end, ringColor, 1f);
            }
        }

        var particles = World.Bodies
            .Where(b => b.UserData is Particle)
            .Select(b => (Particle)b.UserData!)
            .ToList();

        foreach (var particle in particles)
        {
            var relPos = particle.Body.Position - _accelerator.Center;
            double dist = relPos.Magnitude;
             if (dist < _accelerator.Radius * 1.5)
            {
                var relPosVec = new Vector2D(relPos.X, relPos.Y);
                var viewPos = center + relPosVec.Normalized * (radius * 0.5f);

                double energy = particle.GetKineticEnergy();
                float barLength = Math.Min((float)energy * 0.1f, radius * 0.4f);

                var end = viewPos + relPosVec.Normalized * barLength;
                Color4 color = ParticleColors.GetValueOrDefault(particle.Properties.Type, Color4.White);

                DrawLine(viewPos, end, color, 2f);
            }
        }

        DrawCircle(center, radius * 0.05f, Color4.White, true);
    }
}
