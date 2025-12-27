using System;
using System.Collections.Generic;
using System.Diagnostics;
using Artemis;
using Artemis.Core;
using Artemis.Particles;

namespace Artemis.Demo
{
    /// <summary>
    /// Sand Tetris 3D - A physics-based puzzle game using Artemis Physics Engine.
    ///
    /// Gameplay:
    /// - Each turn, you receive a sand ball of random size and color
    /// - Choose where to drop it (X and Z coordinates)
    /// - Sand falls, rolls, and spreads according to physics
    /// - If sand of the same color connects from one side to the other, it disappears
    /// - You lose if the terrarium fills up
    /// - You win if you create a hole at the bottom
    ///
    /// The terrarium is THIN (not a 3D cube) to make it playable for humans.
    /// A full 3D Tetris would be too difficult.
    /// </summary>
    public class SandTetris3D
    {
        #region Constants

        private const double TerrariumWidth = 10.0;
        private const double TerrariumHeight = 15.0;
        private const double TerrariumDepth = 2.5;  // THIN terrarium for playability
        private const double MinBallRadius = 0.6;
        private const double MaxBallRadius = 1.5;
        private const double ParticleRadius = 0.12;
        private const int ScorePerClear = 100;
        private const int BonusPerChain = 50;
        private const float ShineAnimationDuration = 1.5f;  // Seconds before layer disappears

        #endregion

        #region Fields

        private readonly SandSimulation _simulation;
        private readonly Random _random;
        private readonly uint[] _colors;
        private readonly Stopwatch _performanceTimer;
        private int _score;
        private int _turn;
        private int _lastGroupId;
        private bool _gameOver;
        private bool _victory;
        private SandBallInfo _currentBall;

        // Shining layer animation
        private readonly Dictionary<int, float> _shiningParticles;
        private readonly HashSet<int> _particlesToRemove;

        // Performance tracking
        private double _lastUpdateTimeMs;
        private double _averageUpdateTimeMs;
        private int _updateCount;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the current score.
        /// </summary>
        public int Score => _score;

        /// <summary>
        /// Gets the current turn number.
        /// </summary>
        public int Turn => _turn;

        /// <summary>
        /// Gets whether the game is over.
        /// </summary>
        public bool IsGameOver => _gameOver;

        /// <summary>
        /// Gets whether the player has won.
        /// </summary>
        public bool IsVictory => _victory;

        /// <summary>
        /// Gets the current sand ball to drop.
        /// </summary>
        public SandBallInfo CurrentBall => _currentBall;

        /// <summary>
        /// Gets the simulation.
        /// </summary>
        public SandSimulation Simulation => _simulation;

        /// <summary>
        /// Gets the terrarium bounds.
        /// </summary>
        public AABB TerrariumBounds => _simulation.Bounds;

        /// <summary>
        /// Gets the current particle count (real-time grain counter).
        /// </summary>
        public int ParticleCount => _simulation.ParticleCount;

        /// <summary>
        /// Gets the active particle count.
        /// </summary>
        public int ActiveParticleCount => _simulation.ActiveParticleCount;

        /// <summary>
        /// Gets the last update time in milliseconds.
        /// </summary>
        public double LastUpdateTimeMs => _lastUpdateTimeMs;

        /// <summary>
        /// Gets the average update time in milliseconds.
        /// </summary>
        public double AverageUpdateTimeMs => _averageUpdateTimeMs;

        /// <summary>
        /// Gets particles that are currently shining (about to disappear).
        /// Key is particle index, value is shine progress (0-1).
        /// </summary>
        public IReadOnlyDictionary<int, float> ShiningParticles => _shiningParticles;

        /// <summary>
        /// Gets performance stats as a formatted string.
        /// </summary>
        public string PerformanceStats => $"Grains: {ActiveParticleCount} | Update: {_lastUpdateTimeMs:F1}ms | Avg: {_averageUpdateTimeMs:F1}ms";

        /// <summary>
        /// Gets the grain counter dialog text.
        /// </summary>
        public string GrainCounterDialog => $"Sand Grains: {ActiveParticleCount:N0}";

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new Sand Tetris 3D game.
        /// </summary>
        /// <param name="startWithRandomTopography">If true, starts with random disconnected sand layers.</param>
        public SandTetris3D(bool startWithRandomTopography = true)
        {
            var bounds = new AABB(
                new Vector3D(-TerrariumWidth / 2, 0, -TerrariumDepth / 2),
                new Vector3D(TerrariumWidth / 2, TerrariumHeight, TerrariumDepth / 2)
            );

            _simulation = Physics.CreateSandSimulation(bounds, ParticleRadius * 2.5);
            _simulation.ParticleRadius = ParticleRadius;
            _simulation.Friction = 0.35;      // Slightly less friction for better flow
            _simulation.Restitution = 0.08;   // Lower bounce for stability

            _random = new Random();
            _lastGroupId = 0;
            _performanceTimer = new Stopwatch();
            _shiningParticles = new Dictionary<int, float>();
            _particlesToRemove = new HashSet<int>();

            // Define 6 distinct colors (ARGB format)
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
            _updateCount = 0;
            _averageUpdateTimeMs = 0;

            if (startWithRandomTopography)
            {
                GenerateRandomTopography();
            }

            GenerateNextBall();
        }

        #endregion

        #region Random Topography

        /// <summary>
        /// Generates random disconnected sand layers as starting topography.
        /// These create opportunities for connections without making the game too easy.
        /// </summary>
        private void GenerateRandomTopography()
        {
            int numClusters = _random.Next(4, 8);  // 4-7 clusters

            for (int i = 0; i < numClusters; i++)
            {
                uint color = _colors[_random.Next(_colors.Length)];

                // Random position in lower half of terrarium
                double x = (_random.NextDouble() - 0.5) * (TerrariumWidth - 2);
                double y = _random.NextDouble() * (TerrariumHeight * 0.3) + 0.5;  // Lower 30%
                double z = (_random.NextDouble() - 0.5) * (TerrariumDepth - 0.5);

                // Random cluster shape - either horizontal layer or small pile
                bool isHorizontalLayer = _random.NextDouble() > 0.5;

                if (isHorizontalLayer)
                {
                    // Create a thin horizontal layer (not spanning full width - gaps for connections)
                    double layerWidth = TerrariumWidth * (0.2 + _random.NextDouble() * 0.3);  // 20-50% width
                    double layerDepth = TerrariumDepth * 0.8;
                    double layerHeight = ParticleRadius * (3 + _random.Next(4));

                    AddSandLayer(
                        new Vector3D(x, y, z),
                        layerWidth, layerHeight, layerDepth,
                        color
                    );
                }
                else
                {
                    // Create a small pile
                    double pileRadius = 0.4 + _random.NextDouble() * 0.6;
                    _simulation.AddSandBall(
                        new Vector3D(x, y + pileRadius, z),
                        pileRadius,
                        color,
                        ParticleRadius
                    );
                }
            }

            // Let the initial topography settle briefly
            for (int i = 0; i < 10; i++)
            {
                _simulation.Update(0.02);
            }
        }

        /// <summary>
        /// Adds a rectangular layer of sand.
        /// </summary>
        private void AddSandLayer(Vector3D center, double width, double height, double depth, uint color)
        {
            double spacing = ParticleRadius * 2.2;

            int countX = (int)(width / spacing);
            int countY = (int)(height / spacing);
            int countZ = (int)(depth / spacing);

            for (int ix = 0; ix < countX; ix++)
            {
                for (int iy = 0; iy < countY; iy++)
                {
                    for (int iz = 0; iz < countZ; iz++)
                    {
                        double px = center.X - width / 2 + ix * spacing + spacing / 2;
                        double py = center.Y - height / 2 + iy * spacing + spacing / 2;
                        double pz = center.Z - depth / 2 + iz * spacing + spacing / 2;

                        // Add small random offset for natural look
                        px += (_random.NextDouble() - 0.5) * ParticleRadius * 0.3;
                        py += (_random.NextDouble() - 0.5) * ParticleRadius * 0.3;
                        pz += (_random.NextDouble() - 0.5) * ParticleRadius * 0.3;

                        _simulation.AddParticle(new Vector3D(px, py, pz), Vector3D.Zero, color);
                    }
                }
            }
        }

        #endregion

        #region Game Logic

        /// <summary>
        /// Generates the next sand ball to drop.
        /// </summary>
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
                ParticleCount = EstimateParticleCount(radius)
            };
        }

        /// <summary>
        /// Drops the current sand ball at the specified position.
        /// </summary>
        /// <param name="x">X coordinate (-TerrariumWidth/2 to TerrariumWidth/2)</param>
        /// <param name="z">Z coordinate (-TerrariumDepth/2 to TerrariumDepth/2)</param>
        /// <returns>True if the drop was successful.</returns>
        public bool DropBall(double x, double z)
        {
            if (_gameOver)
                return false;

            // Clamp position to within terrarium
            x = Math.Clamp(x, -TerrariumWidth / 2 + _currentBall.Radius,
                          TerrariumWidth / 2 - _currentBall.Radius);
            z = Math.Clamp(z, -TerrariumDepth / 2 + _currentBall.Radius,
                          TerrariumDepth / 2 - _currentBall.Radius);

            // Drop from just below the top
            double y = TerrariumHeight - _currentBall.Radius - 0.5;

            var center = new Vector3D(x, y, z);

            // Create sand ball
            _lastGroupId++;
            var (startIndex, count) = _simulation.AddSandBall(
                center,
                _currentBall.Radius,
                _currentBall.Color,
                ParticleRadius
            );

            return true;
        }

        /// <summary>
        /// Updates the simulation with performance tracking.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void Update(double deltaTime)
        {
            if (_gameOver)
                return;

            _performanceTimer.Restart();

            // Update shining particles animation
            UpdateShiningParticles((float)deltaTime);

            // Run simulation with substeps for stability
            int substeps = deltaTime > 0.02 ? 2 : 1;
            double subDt = deltaTime / substeps;

            for (int i = 0; i < substeps; i++)
            {
                _simulation.Update(subDt);
            }

            _performanceTimer.Stop();
            _lastUpdateTimeMs = _performanceTimer.Elapsed.TotalMilliseconds;

            // Update rolling average
            _updateCount++;
            _averageUpdateTimeMs = _averageUpdateTimeMs + (_lastUpdateTimeMs - _averageUpdateTimeMs) / Math.Min(_updateCount, 100);
        }

        /// <summary>
        /// Updates the shining particle animation.
        /// </summary>
        private void UpdateShiningParticles(float deltaTime)
        {
            _particlesToRemove.Clear();

            foreach (var kvp in _shiningParticles)
            {
                float newProgress = kvp.Value + deltaTime / ShineAnimationDuration;

                if (newProgress >= 1.0f)
                {
                    _particlesToRemove.Add(kvp.Key);
                }
                else
                {
                    _shiningParticles[kvp.Key] = newProgress;
                }
            }

            // Remove particles that finished shining
            foreach (int idx in _particlesToRemove)
            {
                _simulation.RemoveParticle(idx);
                _shiningParticles.Remove(idx);
            }
        }

        /// <summary>
        /// Checks for completed lines and starts the shining animation.
        /// Should be called after simulation settles.
        /// </summary>
        /// <returns>Number of particles that will be cleared.</returns>
        public int CheckAndClearLines()
        {
            int totalToRemove = 0;

            // Check each color
            foreach (uint color in _colors)
            {
                // Check if color spans from left to right (X axis) - primary clear direction
                if (_simulation.CheckSpansAxis(color, 0))
                {
                    var toShine = FindSpanningParticles(color, 0);
                    StartShiningAnimation(toShine);
                    totalToRemove += toShine.Count;
                }
            }

            if (totalToRemove > 0)
            {
                _score += ScorePerClear + (totalToRemove * 10);
            }

            return totalToRemove;
        }

        /// <summary>
        /// Finds all particles of a color that are part of a spanning connection.
        /// </summary>
        private List<int> FindSpanningParticles(uint color, int axis)
        {
            var result = new List<int>();
            var visited = new HashSet<int>();

            for (int i = 0; i < _simulation.Particles.Count; i++)
            {
                if (visited.Contains(i)) continue;

                var p = _simulation.Particles[i];
                if (!p.IsActive || p.Color != color) continue;

                var connected = _simulation.FindConnectedParticles(i);

                foreach (int idx in connected)
                {
                    visited.Add(idx);
                }

                if (CheckGroupSpansAxis(connected, axis))
                {
                    result.AddRange(connected);
                }
            }

            return result;
        }

        /// <summary>
        /// Starts the shining animation for particles about to be removed.
        /// </summary>
        private void StartShiningAnimation(List<int> particleIndices)
        {
            foreach (int idx in particleIndices)
            {
                if (!_shiningParticles.ContainsKey(idx))
                {
                    _shiningParticles[idx] = 0;
                }
            }
        }

        /// <summary>
        /// Gets the shine color for a particle (brightened version of original).
        /// </summary>
        public uint GetShineColor(int particleIndex, uint originalColor)
        {
            if (!_shiningParticles.TryGetValue(particleIndex, out float progress))
            {
                return originalColor;
            }

            // Pulse effect: sine wave oscillation
            float pulse = (float)(Math.Sin(progress * Math.PI * 4) * 0.5 + 0.5);
            float brightness = 0.5f + pulse * 0.5f;

            // Extract RGB
            byte r = (byte)((originalColor >> 16) & 0xFF);
            byte g = (byte)((originalColor >> 8) & 0xFF);
            byte b = (byte)(originalColor & 0xFF);

            // Brighten
            r = (byte)Math.Min(255, r + (int)(255 - r) * brightness);
            g = (byte)Math.Min(255, g + (int)(255 - g) * brightness);
            b = (byte)Math.Min(255, b + (int)(255 - b) * brightness);

            return 0xFF000000 | ((uint)r << 16) | ((uint)g << 8) | b;
        }

        /// <summary>
        /// Checks if a particle is currently shining.
        /// </summary>
        public bool IsParticleShining(int particleIndex)
        {
            return _shiningParticles.ContainsKey(particleIndex);
        }

        private bool CheckGroupSpansAxis(List<int> indices, int axis)
        {
            if (indices.Count == 0) return false;

            double min = double.MaxValue;
            double max = double.MinValue;

            foreach (int idx in indices)
            {
                var p = _simulation.Particles[idx];
                double val = axis switch
                {
                    0 => p.Position.X,
                    1 => p.Position.Y,
                    2 => p.Position.Z,
                    _ => 0
                };
                min = Math.Min(min, val);
                max = Math.Max(max, val);
            }

            double axisSize = axis switch
            {
                0 => TerrariumWidth,
                1 => TerrariumHeight,
                2 => TerrariumDepth,
                _ => 0
            };

            // Require 75% span for a match
            return (max - min) >= axisSize * 0.75;
        }

        /// <summary>
        /// Checks game end conditions.
        /// </summary>
        public void CheckGameState()
        {
            if (_gameOver)
                return;

            // Check win condition: hole at bottom
            if (_simulation.HasHoleAtBottom() && _simulation.ParticleCount > 50)
            {
                _victory = true;
                _gameOver = true;
                return;
            }

            // Check lose condition: terrarium full
            if (_simulation.IsFull())
            {
                _victory = false;
                _gameOver = true;
                return;
            }
        }

        /// <summary>
        /// Advances to the next turn.
        /// </summary>
        public void NextTurn()
        {
            if (!_gameOver)
            {
                GenerateNextBall();
            }
        }

        /// <summary>
        /// Resets the game.
        /// </summary>
        /// <param name="withRandomTopography">If true, generates random starting topography.</param>
        public void Reset(bool withRandomTopography = true)
        {
            _simulation.Clear();
            _shiningParticles.Clear();
            _score = 0;
            _turn = 0;
            _gameOver = false;
            _victory = false;
            _lastGroupId = 0;
            _updateCount = 0;
            _averageUpdateTimeMs = 0;

            if (withRandomTopography)
            {
                GenerateRandomTopography();
            }

            GenerateNextBall();
        }

        #endregion

        #region Helpers

        private int EstimateParticleCount(double ballRadius)
        {
            double volume = (4.0 / 3.0) * Math.PI * ballRadius * ballRadius * ballRadius;
            double particleVolume = (4.0 / 3.0) * Math.PI * ParticleRadius * ParticleRadius * ParticleRadius;
            double packingEfficiency = 0.64; // Random close packing
            return (int)(volume * packingEfficiency / particleVolume);
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

        /// <summary>
        /// Gets the color as RGB components.
        /// </summary>
        public static (byte r, byte g, byte b) GetRGB(uint color)
        {
            return (
                (byte)((color >> 16) & 0xFF),
                (byte)((color >> 8) & 0xFF),
                (byte)(color & 0xFF)
            );
        }

        /// <summary>
        /// Drops the ball at a screen position (for mouse/touch input).
        /// Converts screen coordinates to world coordinates.
        /// </summary>
        /// <param name="screenX">Screen X position (0 to screenWidth).</param>
        /// <param name="screenY">Screen Y position (0 to screenHeight).</param>
        /// <param name="screenWidth">Total screen width.</param>
        /// <param name="screenHeight">Total screen height.</param>
        /// <returns>True if the drop was successful.</returns>
        public bool DropBallAtScreenPosition(double screenX, double screenY, double screenWidth, double screenHeight)
        {
            // Convert screen coordinates to world coordinates
            // Screen (0,0) is top-left, world center is (0,0)
            double normalizedX = screenX / screenWidth;  // 0 to 1
            double normalizedY = screenY / screenHeight; // 0 to 1

            double worldX = (normalizedX - 0.5) * TerrariumWidth;
            double worldZ = (normalizedY - 0.5) * TerrariumDepth;

            return DropBall(worldX, worldZ);
        }

        /// <summary>
        /// Gets the world position from a screen position.
        /// </summary>
        public Vector3D ScreenToWorldPosition(double screenX, double screenY, double screenWidth, double screenHeight, double worldY = 0)
        {
            double normalizedX = screenX / screenWidth;
            double normalizedY = screenY / screenHeight;

            double worldX = (normalizedX - 0.5) * TerrariumWidth;
            double worldZ = (normalizedY - 0.5) * TerrariumDepth;

            return new Vector3D(worldX, worldY, worldZ);
        }

        /// <summary>
        /// Gets the screen position from a world position.
        /// </summary>
        public (double screenX, double screenY) WorldToScreenPosition(Vector3D worldPos, double screenWidth, double screenHeight)
        {
            double normalizedX = (worldPos.X / TerrariumWidth) + 0.5;
            double normalizedY = (worldPos.Z / TerrariumDepth) + 0.5;

            return (normalizedX * screenWidth, normalizedY * screenHeight);
        }

        /// <summary>
        /// Gets the terrarium dimensions for rendering.
        /// </summary>
        public (double width, double height, double depth) GetTerrariumSize()
        {
            return (TerrariumWidth, TerrariumHeight, TerrariumDepth);
        }

        /// <summary>
        /// Gets a formatted string for the real-time grain counter dialog.
        /// </summary>
        public string GetGrainCounterText()
        {
            return $"=== Sand Grains ===\n" +
                   $"Active: {ActiveParticleCount:N0}\n" +
                   $"Total: {ParticleCount:N0}\n" +
                   $"Shining: {_shiningParticles.Count}\n" +
                   $"Update: {_lastUpdateTimeMs:F1}ms";
        }

        #endregion
    }

    /// <summary>
    /// Information about a sand ball to be dropped.
    /// </summary>
    public struct SandBallInfo
    {
        /// <summary>Radius of the sand ball.</summary>
        public double Radius;

        /// <summary>Color of the sand particles (ARGB).</summary>
        public uint Color;

        /// <summary>Name of the color.</summary>
        public string ColorName;

        /// <summary>Estimated number of particles.</summary>
        public int ParticleCount;
    }
}
