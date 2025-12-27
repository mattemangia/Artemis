using System;
using System.Diagnostics;
using System.Threading;
using Artemis.Core;

namespace Artemis.Demo
{
    /// <summary>
    /// Sand Tetris 3D - Console Demo
    /// A physics-based puzzle game demonstrating the Artemis Physics Engine.
    /// </summary>
    class Program
    {
        private static SandTetris3D? _game;
        private static bool _running = true;
        private static double _dropX = 0;
        private static double _dropZ = 0;
        private static bool _waitingForSettle = false;
        private static int _settleFrames = 0;
        private const int SettleFramesRequired = 60;

        static void Main(string[] args)
        {
            Console.Clear();
            Console.CursorVisible = false;

            PrintTitle();

            Console.WriteLine("\n🏖️  SAND TETRIS 3D - Artemis Physics Engine Demo");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            Console.WriteLine("📋 RULES:");
            Console.WriteLine("   • Drop sand balls into the terrarium");
            Console.WriteLine("   • Same-colored sand connecting wall-to-wall disappears");
            Console.WriteLine("   • Create a hole at the bottom to WIN");
            Console.WriteLine("   • Fill up the terrarium and you LOSE\n");

            Console.WriteLine("🎮 CONTROLS:");
            Console.WriteLine("   Arrow Keys: Move drop position (X/Z)");
            Console.WriteLine("   SPACE/ENTER: Drop sand ball");
            Console.WriteLine("   R: Reset game");
            Console.WriteLine("   Q: Quit\n");

            Console.WriteLine("Press any key to start...");
            Console.ReadKey(true);

            StartGame();
        }

        static void PrintTitle()
        {
            Console.ForegroundColor = ConsoleColor.Cyan;
            Console.WriteLine(@"
    ╔═══════════════════════════════════════════════════════════════╗
    ║                                                               ║
    ║     █████╗ ██████╗ ████████╗███████╗███╗   ███╗██╗███████╗   ║
    ║    ██╔══██╗██╔══██╗╚══██╔══╝██╔════╝████╗ ████║██║██╔════╝   ║
    ║    ███████║██████╔╝   ██║   █████╗  ██╔████╔██║██║███████╗   ║
    ║    ██╔══██║██╔══██╗   ██║   ██╔══╝  ██║╚██╔╝██║██║╚════██║   ║
    ║    ██║  ██║██║  ██║   ██║   ███████╗██║ ╚═╝ ██║██║███████║   ║
    ║    ╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝     ╚═╝╚═╝╚══════╝   ║
    ║                                                               ║
    ║              ═══ PHYSICS ENGINE v1.0.0 ═══                    ║
    ║                                                               ║
    ╚═══════════════════════════════════════════════════════════════╝
");
            Console.ResetColor();
        }

        static void StartGame()
        {
            _game = new SandTetris3D();
            _dropX = 0;
            _dropZ = 0;
            _running = true;
            _waitingForSettle = false;

            var stopwatch = Stopwatch.StartNew();
            double lastTime = 0;
            double physicsAccumulator = 0;
            const double PhysicsStep = 1.0 / 60.0;
            const double RenderInterval = 1.0 / 30.0;
            double lastRenderTime = 0;

            while (_running)
            {
                double currentTime = stopwatch.Elapsed.TotalSeconds;
                double deltaTime = currentTime - lastTime;
                lastTime = currentTime;

                // Handle input
                HandleInput();

                // Update physics
                if (_waitingForSettle)
                {
                    physicsAccumulator += deltaTime;
                    while (physicsAccumulator >= PhysicsStep)
                    {
                        _game.Update(PhysicsStep);
                        physicsAccumulator -= PhysicsStep;
                        _settleFrames++;
                    }

                    // Check if settled
                    if (_settleFrames >= SettleFramesRequired)
                    {
                        int cleared = _game.CheckAndClearLines();
                        _game.CheckGameState();

                        if (!_game.IsGameOver)
                        {
                            _game.NextTurn();
                        }

                        _waitingForSettle = false;
                        _settleFrames = 0;
                    }
                }

                // Render
                if (currentTime - lastRenderTime >= RenderInterval)
                {
                    Render();
                    lastRenderTime = currentTime;
                }

                // Check game over
                if (_game.IsGameOver)
                {
                    RenderGameOver();
                    break;
                }

                Thread.Sleep(16); // ~60 fps
            }

            Console.CursorVisible = true;
        }

        static void HandleInput()
        {
            if (!Console.KeyAvailable)
                return;

            var key = Console.ReadKey(true);

            if (_waitingForSettle && key.Key != ConsoleKey.Q && key.Key != ConsoleKey.R)
                return;

            double moveStep = 0.5;

            switch (key.Key)
            {
                case ConsoleKey.LeftArrow:
                    _dropX = Math.Max(_dropX - moveStep, -4.5);
                    break;

                case ConsoleKey.RightArrow:
                    _dropX = Math.Min(_dropX + moveStep, 4.5);
                    break;

                case ConsoleKey.UpArrow:
                    _dropZ = Math.Max(_dropZ - moveStep, -4.5);
                    break;

                case ConsoleKey.DownArrow:
                    _dropZ = Math.Min(_dropZ + moveStep, 4.5);
                    break;

                case ConsoleKey.Spacebar:
                case ConsoleKey.Enter:
                    if (!_waitingForSettle && _game != null && !_game.IsGameOver)
                    {
                        _game.DropBall(_dropX, _dropZ);
                        _waitingForSettle = true;
                        _settleFrames = 0;
                    }
                    break;

                case ConsoleKey.R:
                    _game?.Reset();
                    _dropX = 0;
                    _dropZ = 0;
                    _waitingForSettle = false;
                    break;

                case ConsoleKey.Q:
                    _running = false;
                    break;
            }
        }

        static void Render()
        {
            if (_game == null) return;

            Console.SetCursorPosition(0, 0);

            // Header
            Console.ForegroundColor = ConsoleColor.White;
            Console.WriteLine("╔══════════════════════════════════════════════════════════════════╗");
            Console.WriteLine($"║  SAND TETRIS 3D                   Turn: {_game.Turn,-4} Score: {_game.Score,-8} ║");
            Console.WriteLine("╠══════════════════════════════════════════════════════════════════╣");

            // Current ball info
            var ball = _game.CurrentBall;
            var (r, g, b) = SandTetris3D.GetRGB(ball.Color);
            Console.Write("║  Next Ball: ");
            Console.ForegroundColor = GetConsoleColor(ball.Color);
            Console.Write($"● {ball.ColorName}");
            Console.ForegroundColor = ConsoleColor.White;
            Console.WriteLine($" (Size: {ball.Radius:F1}, ~{ball.ParticleCount} particles)          ║");

            // Drop position
            Console.WriteLine($"║  Drop Position: X={_dropX,5:F1}  Z={_dropZ,5:F1}                              ║");
            Console.WriteLine("╠══════════════════════════════════════════════════════════════════╣");

            // Terrarium view (simplified top-down)
            RenderTerrariumTopDown();

            // Status
            Console.ForegroundColor = ConsoleColor.White;
            Console.WriteLine("╠══════════════════════════════════════════════════════════════════╣");

            if (_waitingForSettle)
            {
                Console.ForegroundColor = ConsoleColor.Yellow;
                Console.WriteLine($"║  ⏳ Sand settling... ({_settleFrames}/{SettleFramesRequired})                               ║");
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine("║  ✓ Ready to drop! Use arrow keys to position, SPACE to drop.    ║");
            }

            Console.ForegroundColor = ConsoleColor.White;
            Console.WriteLine($"║  Particles: {_game.Simulation.ParticleCount,-6}                                          ║");
            Console.WriteLine("╚══════════════════════════════════════════════════════════════════╝");
        }

        static void RenderTerrariumTopDown()
        {
            if (_game == null) return;

            const int GridWidth = 20;
            const int GridHeight = 20;
            char[,] grid = new char[GridHeight, GridWidth];
            ConsoleColor[,] colors = new ConsoleColor[GridHeight, GridWidth];

            // Initialize grid
            for (int z = 0; z < GridHeight; z++)
            {
                for (int x = 0; x < GridWidth; x++)
                {
                    grid[z, x] = '·';
                    colors[z, x] = ConsoleColor.DarkGray;
                }
            }

            // Map particles to grid (top-down view, showing tallest particle at each position)
            var bounds = _game.TerrariumBounds;
            double cellWidth = (bounds.Max.X - bounds.Min.X) / GridWidth;
            double cellDepth = (bounds.Max.Z - bounds.Min.Z) / GridHeight;

            foreach (var particle in _game.Simulation.Particles)
            {
                if (!particle.IsActive) continue;

                int gx = (int)((particle.Position.X - bounds.Min.X) / cellWidth);
                int gz = (int)((particle.Position.Z - bounds.Min.Z) / cellDepth);

                gx = Math.Clamp(gx, 0, GridWidth - 1);
                gz = Math.Clamp(gz, 0, GridHeight - 1);

                // Use density to determine character
                double heightRatio = particle.Position.Y / bounds.Max.Y;
                char c = heightRatio > 0.8 ? '█' : heightRatio > 0.5 ? '▓' : heightRatio > 0.25 ? '▒' : '░';

                grid[gz, gx] = c;
                colors[gz, gx] = GetConsoleColor(particle.Color);
            }

            // Mark drop position
            int dropGx = (int)((_dropX - bounds.Min.X + 5) / cellWidth);
            int dropGz = (int)((_dropZ - bounds.Min.Z + 5) / cellDepth);
            dropGx = Math.Clamp(dropGx, 0, GridWidth - 1);
            dropGz = Math.Clamp(dropGz, 0, GridHeight - 1);

            // Render grid
            for (int z = 0; z < GridHeight; z++)
            {
                Console.ForegroundColor = ConsoleColor.White;
                Console.Write("║  ");

                for (int x = 0; x < GridWidth; x++)
                {
                    if (x == dropGx && z == dropGz && !_waitingForSettle)
                    {
                        Console.ForegroundColor = ConsoleColor.White;
                        Console.Write("◎ ");
                    }
                    else
                    {
                        Console.ForegroundColor = colors[z, x];
                        Console.Write($"{grid[z, x]} ");
                    }
                }

                Console.ForegroundColor = ConsoleColor.White;
                Console.WriteLine("                     ║");
            }
        }

        static void RenderGameOver()
        {
            if (_game == null) return;

            Console.Clear();
            Console.SetCursorPosition(0, 5);

            if (_game.IsVictory)
            {
                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine(@"
    ╔═══════════════════════════════════════════════════╗
    ║                                                   ║
    ║   ██╗   ██╗██╗ ██████╗████████╗ ██████╗ ██████╗   ║
    ║   ██║   ██║██║██╔════╝╚══██╔══╝██╔═══██╗██╔══██╗  ║
    ║   ██║   ██║██║██║        ██║   ██║   ██║██████╔╝  ║
    ║   ╚██╗ ██╔╝██║██║        ██║   ██║   ██║██╔══██╗  ║
    ║    ╚████╔╝ ██║╚██████╗   ██║   ╚██████╔╝██║  ██║  ║
    ║     ╚═══╝  ╚═╝ ╚═════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝  ║
    ║                                                   ║
    ║          🎉 You created a hole! 🎉                ║
    ║                                                   ║
    ╚═══════════════════════════════════════════════════╝
");
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine(@"
    ╔═══════════════════════════════════════════════════╗
    ║                                                   ║
    ║    ██████╗  █████╗ ███╗   ███╗███████╗            ║
    ║   ██╔════╝ ██╔══██╗████╗ ████║██╔════╝            ║
    ║   ██║  ███╗███████║██╔████╔██║█████╗              ║
    ║   ██║   ██║██╔══██║██║╚██╔╝██║██╔══╝              ║
    ║   ╚██████╔╝██║  ██║██║ ╚═╝ ██║███████╗            ║
    ║    ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝            ║
    ║                                                   ║
    ║    ██████╗ ██╗   ██╗███████╗██████╗               ║
    ║   ██╔═══██╗██║   ██║██╔════╝██╔══██╗              ║
    ║   ██║   ██║██║   ██║█████╗  ██████╔╝              ║
    ║   ██║   ██║╚██╗ ██╔╝██╔══╝  ██╔══██╗              ║
    ║   ╚██████╔╝ ╚████╔╝ ███████╗██║  ██║              ║
    ║    ╚═════╝   ╚═══╝  ╚══════╝╚═╝  ╚═╝              ║
    ║                                                   ║
    ║          😞 The terrarium is full! 😞             ║
    ║                                                   ║
    ╚═══════════════════════════════════════════════════╝
");
            }

            Console.ForegroundColor = ConsoleColor.Yellow;
            Console.WriteLine($"\n    Final Score: {_game.Score}");
            Console.WriteLine($"    Turns Played: {_game.Turn}");
            Console.WriteLine($"    Particles: {_game.Simulation.ParticleCount}");

            Console.ForegroundColor = ConsoleColor.White;
            Console.WriteLine("\n    Press R to restart or Q to quit...");

            while (true)
            {
                var key = Console.ReadKey(true);
                if (key.Key == ConsoleKey.R)
                {
                    StartGame();
                    break;
                }
                if (key.Key == ConsoleKey.Q)
                {
                    break;
                }
            }
        }

        static ConsoleColor GetConsoleColor(uint argb)
        {
            return argb switch
            {
                0xFFFF4444 => ConsoleColor.Red,
                0xFF44FF44 => ConsoleColor.Green,
                0xFF4444FF => ConsoleColor.Blue,
                0xFFFFFF44 => ConsoleColor.Yellow,
                0xFFFF44FF => ConsoleColor.Magenta,
                0xFF44FFFF => ConsoleColor.Cyan,
                _ => ConsoleColor.White
            };
        }
    }
}
