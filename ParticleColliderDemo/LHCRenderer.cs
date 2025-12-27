using ArtemisEngine;

namespace ParticleColliderDemo;

/// <summary>
/// Renderer for particle collider visualization
/// Inspired by LHC detector event displays (ATLAS, CMS)
/// </summary>
public class LHCRenderer
{
    private int _width;
    private int _height;
    private char[,] _buffer;
    private ConsoleColor[,] _colorBuffer;

    // View parameters
    private Vector2 _cameraCenter;
    private float _zoom;

    public LHCRenderer(int width, int height)
    {
        _width = width;
        _height = height;
        _buffer = new char[width, height];
        _colorBuffer = new ConsoleColor[width, height];
        _cameraCenter = new Vector2(0, 0);
        _zoom = 1.0f;
    }

    public void SetCamera(Vector2 center, float zoom)
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
            float angle1 = (float)i / segments * MathF.PI * 2;
            float angle2 = (float)(i + 1) / segments * MathF.PI * 2;

            Vector2 p1 = accelerator.Center + new Vector2(
                MathF.Cos(angle1) * accelerator.Radius,
                MathF.Sin(angle1) * accelerator.Radius
            );
            Vector2 p2 = accelerator.Center + new Vector2(
                MathF.Cos(angle2) * accelerator.Radius,
                MathF.Sin(angle2) * accelerator.Radius
            );

            DrawLine(p1, p2, '═', ConsoleColor.DarkGray);
        }

        // Draw interaction points
        var ips = accelerator.GetInteractionPoints();
        foreach (var ip in ips)
        {
            DrawCircle(ip, 0.5f, '◉', ConsoleColor.Red);
        }
    }

    public void DrawDetectors(List<Detector> detectors)
    {
        foreach (var detector in detectors)
        {
            char symbol = detector.Type switch
            {
                DetectorType.Tracker => '░',
                DetectorType.Calorimeter => '▒',
                DetectorType.MuonChamber => '▓',
                DetectorType.InteractionPoint => '◉',
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
                ConsoleColor trailColor = particle.Properties.Color;
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
                _ => '?'
            };

            _buffer[screenPos.x, screenPos.y] = symbol;
            _colorBuffer[screenPos.x, screenPos.y] = particle.Properties.Color;
        }

        // Draw velocity vector
        if (particle.Body.Velocity.Length > 1)
        {
            Vector2 end = particle.Body.Position + particle.Body.Velocity.Normalized * 2;
            DrawLine(particle.Body.Position, end, '→', particle.Properties.Color);
        }
    }

    public void DrawCollisionEvent(ParticleCollisionEvent evt, float age)
    {
        if (age > 2.0f) return; // Fade out old events

        // Draw explosion effect at collision point
        float radius = age * 3;
        ConsoleColor color = age < 0.5f ? ConsoleColor.White : ConsoleColor.Yellow;

        DrawCircle(evt.CollisionPoint, radius, '*', color);

        // Draw energy text
        var screenPos = WorldToScreen(evt.CollisionPoint);
        if (IsOnScreen(screenPos) && screenPos.y > 0)
        {
            string energyText = $"{evt.CenterOfMassEnergy:F1}GeV";
            DrawText(screenPos.x, screenPos.y - 1, energyText, ConsoleColor.Yellow);
        }
    }

    private void DrawCircle(Vector2 center, float radius, char symbol, ConsoleColor color)
    {
        int segments = (int)(radius * 20);
        segments = Math.Clamp(segments, 8, 100);

        for (int i = 0; i < segments; i++)
        {
            float angle = (float)i / segments * MathF.PI * 2;
            Vector2 point = center + new Vector2(
                MathF.Cos(angle) * radius,
                MathF.Sin(angle) * radius
            );

            var screenPos = WorldToScreen(point);
            if (IsOnScreen(screenPos))
            {
                _buffer[screenPos.x, screenPos.y] = symbol;
                _colorBuffer[screenPos.x, screenPos.y] = color;
            }
        }
    }

    private void DrawLine(Vector2 start, Vector2 end, char symbol, ConsoleColor color)
    {
        Vector2 delta = end - start;
        float distance = delta.Length;
        int steps = (int)(distance * 2);

        if (steps == 0) return;

        for (int i = 0; i <= steps; i++)
        {
            float t = (float)i / steps;
            Vector2 point = start + delta * t;

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

    private (int x, int y) WorldToScreen(Vector2 worldPos)
    {
        Vector2 offset = worldPos - _cameraCenter;
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
            // Don't use WriteLine - use SetCursorPosition to avoid scrolling
        }

        Console.ResetColor();
    }

    public void DrawUI(CollisionExperiment experiment, int particleCount)
    {
        // Draw compact UI status bar at fixed position (avoid scrolling)
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
        foreach (var type in types)
        {
            var props = ParticleProperties.Create(type);
            Console.ForegroundColor = props.Color;
            Console.Write($"{props.Symbol} ");
            Console.ForegroundColor = ConsoleColor.Gray;
            Console.Write($"= {type,-10} ");

            if ((int)type % 3 == 2)
                Console.WriteLine();
        }

        Console.ResetColor();
        Console.WriteLine();
    }
}
