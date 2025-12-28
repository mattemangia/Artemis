using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Artemis;
using Artemis.Core;
using Artemis.Physics2D;
using Artemis.Materials;

namespace Artemis.Demo
{
    /// <summary>
    /// Sand Tetris 3D - A physics-based puzzle game using Artemis Physics Engine (2D Physics Logic).
    ///
    /// The user requested "parallelepipeds or small cubes with the depth of the sandbox".
    /// This is implemented using a 2D PhysicsWorld but rendering the 2D shapes as 3D blocks
    /// spanning the full depth of the terrarium.
    /// </summary>
    public class SandTetris3D
    {
        #region Constants

        private const double TerrariumWidth = 10.0;
        private const double TerrariumHeight = 15.0;
        private const double TerrariumDepth = 2.5;  // Rendered depth
        private const double MinBallRadius = 0.4;
        private const double MaxBallRadius = 1.2;
        private const double BlockSize = 0.35; // Size of the cubes/parallelepipeds
        private const int ScorePerClear = 100;
        private const int BonusPerChain = 50;
        private const float ShineAnimationDuration = 1.5f;

        #endregion

        #region Fields

        private readonly PhysicsWorld2D _world;
        private readonly Random _random;
        private readonly uint[] _colors;
        private readonly Stopwatch _performanceTimer;
        private int _score;
        private int _turn;
        private bool _gameOver;
        private bool _victory;
        private SandBallInfo _currentBall;

        // Shining animation for cleared blocks
        private readonly Dictionary<string, float> _shiningBodies;
        private readonly HashSet<string> _bodiesToRemove;

        // Performance tracking
        private double _lastUpdateTimeMs;
        private double _averageUpdateTimeMs;
        private int _updateCount;

        #endregion

        #region Properties

        public int Score => _score;
        public int Turn => _turn;
        public bool IsGameOver => _gameOver;
        public bool IsVictory => _victory;
        public SandBallInfo CurrentBall => _currentBall;

        // Expose the PhysicsWorld2D for rendering
        public PhysicsWorld2D World => _world;

        public int ActiveParticleCount => _world.Bodies.Count(b => b.BodyType == BodyType2D.Dynamic);

        public double LastUpdateTimeMs => _lastUpdateTimeMs;
        public double AverageUpdateTimeMs => _averageUpdateTimeMs;
        public string PerformanceStats => $"Blocks: {ActiveParticleCount} | Update: {_lastUpdateTimeMs:F1}ms | Avg: {_averageUpdateTimeMs:F1}ms";
        public string GrainCounterDialog => $"Blocks: {ActiveParticleCount:N0}";

        // Helper to expose shining bodies to renderer
        public IReadOnlyDictionary<string, float> ShiningBodies => _shiningBodies;

        #endregion

        #region Constructor

        public SandTetris3D(bool startWithRandomTopography = true)
        {
            // Initialize 2D world
            _world = Physics.CreateWorld2D();
            _world.Gravity = new Artemis.Physics2D.Vector2D(0, -9.81);
            _world.VelocityIterations = 6;
            _world.PositionIterations = 3;

            // Create container walls
            CreateContainer();

            _random = new Random();
            _performanceTimer = new Stopwatch();
            _shiningBodies = new Dictionary<string, float>();
            _bodiesToRemove = new HashSet<string>();

            // Define colors
            _colors = new uint[]
            {
                0xFFFF4444, // Red
                0xFF44FF44, // Green
                0xFF4444FF, // Blue
                0xFFFFFF44, // Yellow
                0xFFFF44FF, // Magenta
                0xFF44FFFF  // Cyan
            };

            _score = 0;
            _turn = 0;
            _gameOver = false;
            _victory = false;

            if (startWithRandomTopography)
            {
                GenerateRandomTopography();
            }

            GenerateNextBall();
        }

        private void CreateContainer()
        {
            // Floor
            var floor = Physics.CreateStaticBox2D(
                new Artemis.Physics2D.Vector2D(0, -0.5),
                TerrariumWidth,
                1.0);
            floor.Friction = 0.5;
            _world.AddBody(floor);

            // Left Wall
            var leftWall = Physics.CreateStaticBox2D(
                new Artemis.Physics2D.Vector2D(-TerrariumWidth / 2 - 0.5, TerrariumHeight / 2),
                1.0,
                TerrariumHeight + 2.0);
            leftWall.Friction = 0.1;
            _world.AddBody(leftWall);

            // Right Wall
            var rightWall = Physics.CreateStaticBox2D(
                new Artemis.Physics2D.Vector2D(TerrariumWidth / 2 + 0.5, TerrariumHeight / 2),
                1.0,
                TerrariumHeight + 2.0);
            rightWall.Friction = 0.1;
            _world.AddBody(rightWall);
        }

        #endregion

        #region Random Topography

        private void GenerateRandomTopography()
        {
            // Fill bottom 40% with random blocks
            double fillHeight = TerrariumHeight * 0.4;
            int rows = (int)(fillHeight / BlockSize);
            int cols = (int)(TerrariumWidth / BlockSize);

            double startX = -TerrariumWidth / 2 + BlockSize / 2;
            double startY = BlockSize / 2;

        // Use Perlin-like coherent noise for colors to create clusters
        int colorSeed = _random.Next(1000);

            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    // Sine wave pattern
                    double x = startX + c * BlockSize;
                    double waveHeight = Math.Sin(x * 0.5) * 2.0 + 2.0;

                    if (r * BlockSize < waveHeight)
                    {
                    // Ensure the bottom few rows are solid to prevent immediate "hole" victory
                    // and give the player a foundation to clear.
                    bool isSolidLayer = r < 3;

                    // Add some jitter for upper layers
                    if (isSolidLayer || _random.NextDouble() > 0.2)
                        {
                        // Use coherent coloring: simple pattern based on position
                        // Modifying the pattern to create larger blobs/layers
                        int colorIndex = (int)((Math.Sin(c * 0.3 + colorSeed) + Math.Cos(r * 0.3)) * 2 + 3) % _colors.Length;
                        if (colorIndex < 0) colorIndex += _colors.Length;

                        uint color = _colors[colorIndex];

                            var body = CreateBlock(x, startY + r * BlockSize, color);
                            // Settle them initially
                            body.IsSleeping = true;
                        }
                    }
                }
            }
        }

        #endregion

        #region Game Logic

        private void GenerateNextBall()
        {
            _turn++;
            double radius = MinBallRadius + _random.NextDouble() * (MaxBallRadius - MinBallRadius);
            uint color = _colors[_random.Next(_colors.Length)];

            _currentBall = new SandBallInfo
            {
                Radius = radius,
                Color = color,
                ColorName = GetColorName(color),
                ParticleCount = (int)((Math.PI * radius * radius) / (BlockSize * BlockSize))
            };
        }

        public bool DropBall(double x, double z)
        {
            if (_gameOver) return false;

            // Clamp X
            x = Math.Clamp(x, -TerrariumWidth / 2 + _currentBall.Radius, TerrariumWidth / 2 - _currentBall.Radius);

            // Start high up
            double y = TerrariumHeight - _currentBall.Radius - 1.0;

            // Spawn a cluster of blocks approximating the circle
            double r = _currentBall.Radius;
            int count = 0;

            // Grid-based packing for the "ball"
            for (double bx = -r; bx <= r; bx += BlockSize)
            {
                for (double by = -r; by <= r; by += BlockSize)
                {
                    if (bx * bx + by * by <= r * r)
                    {
                        CreateBlock(x + bx, y + by, _currentBall.Color);
                        count++;
                    }
                }
            }

            return true;
        }

        private RigidBody2D CreateBlock(double x, double y, uint color)
        {
            // Add some random offset to prevent perfect stacking issues
            double jitter = BlockSize * 0.05;
            x += (_random.NextDouble() - 0.5) * jitter;
            y += (_random.NextDouble() - 0.5) * jitter;

            var body = Physics.CreateBox2D(new Artemis.Physics2D.Vector2D(x, y), BlockSize, BlockSize);
            body.Restitution = 0.1;
            body.Friction = 0.4;
            body.UserData = color; // Store color in UserData
            _world.AddBody(body);
            return body;
        }

        public void Update(double deltaTime)
        {
            if (_gameOver) return;

            _performanceTimer.Restart();

            // Update shining animation
            UpdateShiningBodies((float)deltaTime);

            // Step physics
            _world.Step(deltaTime);

            _performanceTimer.Stop();
            _lastUpdateTimeMs = _performanceTimer.Elapsed.TotalMilliseconds;

            _updateCount++;
            _averageUpdateTimeMs = _averageUpdateTimeMs + (_lastUpdateTimeMs - _averageUpdateTimeMs) / Math.Min(_updateCount, 100);
        }

        private void UpdateShiningBodies(float deltaTime)
        {
            _bodiesToRemove.Clear();
            var keys = _shiningBodies.Keys.ToList();

            foreach (var key in keys)
            {
                float newProgress = _shiningBodies[key] + deltaTime / ShineAnimationDuration;
                if (newProgress >= 1.0f)
                {
                    _bodiesToRemove.Add(key);
                }
                else
                {
                    _shiningBodies[key] = newProgress;
                }
            }

            foreach (var id in _bodiesToRemove)
            {
                var body = _world.GetBody(id);
                if (body != null)
                {
                    _world.RemoveBody(body);
                }
                _shiningBodies.Remove(id);
            }
        }

        public int CheckAndClearLines()
        {
            // For blocks, we check if connected blocks of same color span left-to-right
            int totalRemoved = 0;

            foreach (uint color in _colors)
            {
                if (CheckSpansAxis(color))
                {
                    var connected = FindSpanningBodies(color);
                    StartShiningAnimation(connected);
                    totalRemoved += connected.Count;
                }
            }

            if (totalRemoved > 0)
            {
                _score += ScorePerClear + (totalRemoved * 10);
            }

            return totalRemoved;
        }

        private bool CheckSpansAxis(uint color)
        {
            // Find left-most blocks of this color
            double minX = -TerrariumWidth / 2 + BlockSize;
            double maxX = TerrariumWidth / 2 - BlockSize;

            var bodies = _world.Bodies.Where(b => b.BodyType == BodyType2D.Dynamic && b.IsActive && b.UserData is uint c && c == color).ToList();

            var startNodes = bodies.Where(b => b.Position.X <= minX + BlockSize).ToList();
            if (startNodes.Count == 0) return false;

            // BFS to see if we can reach maxX
            var visited = new HashSet<RigidBody2D>();
            var queue = new Queue<RigidBody2D>();

            foreach (var node in startNodes)
            {
                queue.Enqueue(node);
                visited.Add(node);
            }

            while (queue.Count > 0)
            {
                var current = queue.Dequeue();
                if (current.Position.X >= maxX - BlockSize) return true;

                // Find neighbors
                foreach (var other in bodies)
                {
                    if (visited.Contains(other)) continue;

                    double distSq = (current.Position - other.Position).MagnitudeSquared;
                    // BlockSize is diagonal approx ~1.41 * size.
                    // Allow slight gap
                    double maxDist = BlockSize * 1.5;
                    if (distSq < maxDist * maxDist)
                    {
                        visited.Add(other);
                        queue.Enqueue(other);
                    }
                }
            }

            return false;
        }

        private List<RigidBody2D> FindSpanningBodies(uint color)
        {
            // Similar to CheckSpansAxis but returns the connected component that spans
            double minX = -TerrariumWidth / 2 + BlockSize;
            double maxX = TerrariumWidth / 2 - BlockSize;

            var bodies = _world.Bodies.Where(b => b.BodyType == BodyType2D.Dynamic && b.IsActive && b.UserData is uint c && c == color).ToList();
            var startNodes = bodies.Where(b => b.Position.X <= minX + BlockSize).ToList();

            var spanningSet = new HashSet<RigidBody2D>();

            foreach (var startNode in startNodes)
            {
                var component = new HashSet<RigidBody2D>();
                var queue = new Queue<RigidBody2D>();
                queue.Enqueue(startNode);
                component.Add(startNode);

                bool reachesEnd = false;

                while (queue.Count > 0)
                {
                    var current = queue.Dequeue();
                    if (current.Position.X >= maxX - BlockSize) reachesEnd = true;

                    foreach (var other in bodies)
                    {
                        if (component.Contains(other)) continue;

                        double distSq = (current.Position - other.Position).MagnitudeSquared;
                        double maxDist = BlockSize * 1.5;

                        if (distSq < maxDist * maxDist)
                        {
                            component.Add(other);
                            queue.Enqueue(other);
                        }
                    }
                }

                if (reachesEnd)
                {
                    foreach (var b in component) spanningSet.Add(b);
                }
            }

            return spanningSet.ToList();
        }

        private void StartShiningAnimation(List<RigidBody2D> bodies)
        {
            foreach (var body in bodies)
            {
                if (!_shiningBodies.ContainsKey(body.Id))
                {
                    _shiningBodies[body.Id] = 0;
                }
            }
        }

        public void CheckGameState()
        {
            if (_gameOver) return;

            // Lose if full
            // Check if any settled body is above a certain height
            double topLimit = TerrariumHeight - 2.0;
            if (_world.Bodies.Any(b => b.BodyType == BodyType2D.Dynamic && b.IsActive && b.Position.Y > topLimit && Math.Abs(b.Velocity.Y) < 0.1))
            {
                _gameOver = true;
                _victory = false;
                return;
            }

            // Win if hole at bottom (and enough blocks)
            // Ensure we have played a few turns to avoid instant win on start glitches
            if (ActiveParticleCount > 50 && _turn > 5)
            {
                // Check bottom row for gaps
                // Simple check: Raycast or interval check across the bottom
                bool hasGap = false;
                double yCheck = 0.5; // Just above floor

                // Sweep X
                for (double x = -TerrariumWidth/2 + 1.0; x < TerrariumWidth/2 - 1.0; x += BlockSize)
                {
                    var bodiesAtX = _world.QueryPoint(new Artemis.Physics2D.Vector2D(x, yCheck));
                    if (bodiesAtX.Count == 0)
                    {
                        hasGap = true;
                        break;
                    }
                }

                if (hasGap)
                {
                    _victory = true;
                    _gameOver = true;
                }
            }
        }

        public void NextTurn()
        {
            if (!_gameOver) GenerateNextBall();
        }

        public void Reset(bool withRandomTopography = true)
        {
            _world.ClearBodies();
            CreateContainer(); // Re-add walls

            _shiningBodies.Clear();
            _score = 0;
            _turn = 0;
            _gameOver = false;
            _victory = false;

            if (withRandomTopography)
            {
                GenerateRandomTopography();
            }

            GenerateNextBall();
        }

        #endregion

        #region Helpers

        public (double width, double height, double depth) GetTerrariumSize()
        {
            return (TerrariumWidth, TerrariumHeight, TerrariumDepth);
        }

        public uint GetShineColor(uint originalColor, float progress)
        {
            float pulse = (float)(Math.Sin(progress * Math.PI * 4) * 0.5 + 0.5);
            float brightness = 0.5f + pulse * 0.5f;

            byte r = (byte)((originalColor >> 16) & 0xFF);
            byte g = (byte)((originalColor >> 8) & 0xFF);
            byte b = (byte)(originalColor & 0xFF);

            r = (byte)Math.Min(255, r + (int)(255 - r) * brightness);
            g = (byte)Math.Min(255, g + (int)(255 - g) * brightness);
            b = (byte)Math.Min(255, b + (int)(255 - b) * brightness);

            return 0xFF000000 | ((uint)r << 16) | ((uint)g << 8) | b;
        }

        private string GetColorName(uint color)
        {
            return color switch
            {
                0xFFFF4444 => "Red",
                0xFF44FF44 => "Green",
                0xFF4444FF => "Blue",
                0xFFFFFF44 => "Yellow",
                0xFFFF44FF => "Magenta",
                0xFF44FFFF => "Cyan",
                _ => "Unknown"
            };
        }

        public static (byte r, byte g, byte b) GetRGB(uint color)
        {
            return (
                (byte)((color >> 16) & 0xFF),
                (byte)((color >> 8) & 0xFF),
                (byte)(color & 0xFF)
            );
        }

        #endregion
    }

    public struct SandBallInfo
    {
        public double Radius;
        public uint Color;
        public string ColorName;
        public int ParticleCount;
    }
}
