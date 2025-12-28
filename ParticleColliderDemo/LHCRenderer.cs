using Artemis.Physics2D;
using Artemis.Physics2D.Particles;

namespace ParticleColliderDemo;

/// <summary>
/// Renderer for particle collider visualization.
/// Inspired by LHC detector event displays (ATLAS, CMS).
/// </summary>
public class LHCRenderer
{
    private int _width;
    private int _height;
    private char[,] _buffer;
    private ConsoleColor[,] _colorBuffer;

    // View parameters
    private Vector2D _cameraCenter;
    private double _zoom;

    public LHCRenderer(int width, int height)
    {
        _width = width;
        _height = height;
        _buffer = new char[width, height];
        _colorBuffer = new ConsoleColor[width, height];
        _cameraCenter = Vector2D.Zero;
        _zoom = 1.0;
    }

    public void SetCamera(Vector2D center, double zoom)
    {
        _cameraCenter = center;
        _zoom = zoom;
    }

    public void Clear()
    {
        for (int y = 0; y < _height; y++)
        {
            for (int x = 0; x < _width; x++)
            {
                _buffer[x, y] = ' ';
                _colorBuffer[x, y] = ConsoleColor.Black;
            }
        }
    }

    public void DrawAcceleratorRing(Accelerator accelerator)
    {
        // Draw the ring
        int segments = 120;
        for (int i = 0; i < segments; i++)
        {
            double angle1 = (double)i / segments * Math.PI * 2;
            double angle2 = (double)(i + 1) / segments * Math.PI * 2;

            Vector2D p1 = accelerator.Center + Vector2D.FromAngle(angle1, accelerator.Radius);
            Vector2D p2 = accelerator.Center + Vector2D.FromAngle(angle2, accelerator.Radius);

            DrawLine(p1, p2, '═', ConsoleColor.DarkGray);
        }

        // Draw interaction points
        var ips = accelerator.GetInteractionPoints();
        foreach (var ip in ips)
        {
            DrawCircle(ip, 0.5, '◉', ConsoleColor.Red);
        }
    }

    public void DrawDetectors(List<Detector> detectors)
    {
        foreach (var detector in detectors)
        {
            char symbol = detector.Type switch
            {
                DetectorType2D.Tracker => '░',
                DetectorType2D.Calorimeter => '▒',
                DetectorType2D.MuonChamber => '▓',
                DetectorType2D.InteractionPoint => '◉',
                _ => '?'
            };

            DrawCircle(detector.Position, detector.Radius, symbol, detector.Color);
        }
    }

    public void DrawParticle(Particle particle)
    {
        // Draw trail first
        if (particle.Trail.Count > 1)
        {
            for (int i = 0; i < particle.Trail.Count - 1; i++)
            {
                // Fade trail
                ConsoleColor trailColor = particle.Properties.GetConsoleColor();
                if (i < particle.Trail.Count / 2)
                {
                    trailColor = ConsoleColor.DarkGray;
                }

                DrawLine(particle.Trail[i], particle.Trail[i + 1], '·', trailColor);
            }
        }

        // Draw particle symbol
        var screenPos = WorldToScreen(particle.Body.Position);
        if (IsOnScreen(screenPos))
        {
            char symbol = particle.Properties.Type switch
            {
                ParticleType.Proton => 'P',
                ParticleType.Electron => 'e',
                ParticleType.Positron => '+',
                ParticleType.Neutron => 'N',
                ParticleType.Photon => '☼',
                ParticleType.Higgs => 'H',
                ParticleType.Quark => 'q',
                ParticleType.Muon => 'μ',
                ParticleType.Neutrino => 'ν',
                ParticleType.Pion => 'π',
                _ => '?'
            };

            _buffer[screenPos.x, screenPos.y] = symbol;
            _colorBuffer[screenPos.x, screenPos.y] = particle.Properties.GetConsoleColor();
        }

        // Draw velocity vector
        if (particle.Body.Velocity.Magnitude > 1)
        {
            Vector2D end = particle.Body.Position + particle.Body.Velocity.Normalized * 2;
            DrawLine(particle.Body.Position, end, '→', particle.Properties.GetConsoleColor());
        }
    }

    public void DrawCollisionEvent(ParticleCollisionEvent evt, double age)
    {
        if (age > 2.0) return; // Fade out old events

        // Draw explosion effect at collision point
        double radius = age * 3;
        ConsoleColor color = age < 0.5 ? ConsoleColor.White : ConsoleColor.Yellow;

        DrawCircle(evt.CollisionPoint, radius, '*', color);

        // Draw energy text
        var screenPos = WorldToScreen(evt.CollisionPoint);
        if (IsOnScreen(screenPos) && screenPos.y > 0)
        {
            string energyText = $"{evt.CenterOfMassEnergy:F1}GeV";
            DrawText(screenPos.x, screenPos.y - 1, energyText, ConsoleColor.Yellow);
        }
    }

    private void DrawCircle(Vector2D center, double radius, char symbol, ConsoleColor color)
    {
        int segments = (int)(radius * 20);
        segments = Math.Clamp(segments, 8, 100);

        for (int i = 0; i < segments; i++)
        {
            double angle = (double)i / segments * Math.PI * 2;
            Vector2D point = center + Vector2D.FromAngle(angle, radius);

            var screenPos = WorldToScreen(point);
            if (IsOnScreen(screenPos))
            {
                _buffer[screenPos.x, screenPos.y] = symbol;
                _colorBuffer[screenPos.x, screenPos.y] = color;
            }
        }
    }

    private void DrawLine(Vector2D start, Vector2D end, char symbol, ConsoleColor color)
    {
        Vector2D delta = end - start;
        double distance = delta.Magnitude;
        int steps = (int)(distance * 2);

        if (steps == 0) return;

        for (int i = 0; i <= steps; i++)
        {
            double t = (double)i / steps;
            Vector2D point = start + delta * t;

            var screenPos = WorldToScreen(point);
            if (IsOnScreen(screenPos))
            {
                _buffer[screenPos.x, screenPos.y] = symbol;
                _colorBuffer[screenPos.x, screenPos.y] = color;
            }
        }
    }

    private void DrawText(int x, int y, string text, ConsoleColor color)
    {
        for (int i = 0; i < text.Length; i++)
        {
            int screenX = x + i;
            if (screenX >= 0 && screenX < _width && y >= 0 && y < _height)
            {
                _buffer[screenX, y] = text[i];
                _colorBuffer[screenX, y] = color;
            }
        }
    }

    private (int x, int y) WorldToScreen(Vector2D worldPos)
    {
        Vector2D offset = worldPos - _cameraCenter;
        int x = (int)((offset.X * _zoom) + _width / 2);
        int y = (int)((-offset.Y * _zoom) + _height / 2); // Flip Y

        return (x, y);
    }

    private bool IsOnScreen((int x, int y) screenPos)
    {
        return screenPos.x >= 0 && screenPos.x < _width &&
               screenPos.y >= 0 && screenPos.y < _height;
    }

    public void Render()
    {
        Console.SetCursorPosition(0, 0);

        for (int y = 0; y < _height; y++)
        {
            Console.SetCursorPosition(0, y);
            for (int x = 0; x < _width; x++)
            {
                Console.ForegroundColor = _colorBuffer[x, y];
                Console.Write(_buffer[x, y]);
            }
        }

        Console.ResetColor();
    }

    public void DrawUI(CollisionExperiment experiment, int particleCount)
    {
        // Draw compact UI status bar at fixed position
        int uiLine = _height;

        Console.SetCursorPosition(0, uiLine);
        Console.ForegroundColor = ConsoleColor.Cyan;

        var (b1, b2, avgE) = experiment.Accelerator.GetStatistics();
        string status = $"B1:{b1,3} B2:{b2,3} | E:{avgE:F0}GeV | Col:{experiment.TotalCollisions,4} | Parts:{particleCount,3}";
        Console.Write(status.PadRight(_width));

        Console.SetCursorPosition(0, uiLine + 1);
        Console.ForegroundColor = ConsoleColor.Yellow;
        Console.Write("[SPACE]Inject [1-3]Types [C]Clear [E]Export [Q]Quit".PadRight(_width));

        Console.ResetColor();
    }

    public void DrawParticleKey()
    {
        Console.SetCursorPosition(0, _height + 15);
        Console.WriteLine("Particle Types:");

        var types = Enum.GetValues<ParticleType>();
        int count = 0;
        foreach (var type in types)
        {
            var props = ParticleProperties.Create(type);
            Console.ForegroundColor = props.GetConsoleColor();
            Console.Write($"{props.Symbol} ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.Write($"= {type,-10} ");

            if (++count % 3 == 0)
                Console.WriteLine();
        }

        Console.ResetColor();
        Console.WriteLine();
    }
}
