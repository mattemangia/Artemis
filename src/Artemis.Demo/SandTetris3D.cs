using System;
using System.Collections.Generic;
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
    /// </summary>
    public class SandTetris3D
    {
        #region Constants

        private const double TerrariumWidth = 10.0;
        private const double TerrariumHeight = 15.0;
        private const double TerrariumDepth = 10.0;
        private const double MinBallRadius = 0.8;
        private const double MaxBallRadius = 2.0;
        private const double ParticleRadius = 0.15;
        private const int ScorePerClear = 100;
        private const int BonusPerChain = 50;

        #endregion

        #region Fields

        private readonly SandSimulation _simulation;
        private readonly Random _random;
        private readonly uint[] _colors;
        private int _score;
        private int _turn;
        private int _lastGroupId;
        private bool _gameOver;
        private bool _victory;
        private SandBallInfo _currentBall;

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

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new Sand Tetris 3D game.
        /// </summary>
        public SandTetris3D()
        {
            var bounds = new AABB(
                new Vector3D(-TerrariumWidth / 2, 0, -TerrariumDepth / 2),
                new Vector3D(TerrariumWidth / 2, TerrariumHeight, TerrariumDepth / 2)
            );

            _simulation = Physics.CreateSandSimulation(bounds, ParticleRadius * 3);
            _simulation.ParticleRadius = ParticleRadius;
            _simulation.Friction = 0.4;
            _simulation.Restitution = 0.1;

            _random = new Random();
            _lastGroupId = 0;

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

            GenerateNextBall();
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

            // Assign group ID to all new particles
            for (int i = startIndex; i < startIndex + count && i < _simulation.Particles.Count; i++)
            {
                var p = _simulation.Particles[i];
                // Note: Since SandParticle is a struct, we need to handle this differently
                // For now, we track the last dropped particles by index range
            }

            return true;
        }

        /// <summary>
        /// Updates the simulation.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void Update(double deltaTime)
        {
            if (_gameOver)
                return;

            // Run simulation
            _simulation.Update(deltaTime);
        }

        /// <summary>
        /// Checks for completed lines and clears them.
        /// Should be called after simulation settles.
        /// </summary>
        /// <returns>Number of particles cleared.</returns>
        public int CheckAndClearLines()
        {
            int totalCleared = 0;

            // Check each color
            foreach (uint color in _colors)
            {
                // Check if color spans from left to right (X axis)
                if (_simulation.CheckSpansAxis(color, 0))
                {
                    // Find and remove all connected particles of this color
                    var toRemove = new List<int>();

                    for (int i = 0; i < _simulation.Particles.Count; i++)
                    {
                        var p = _simulation.Particles[i];
                        if (p.IsActive && p.Color == color)
                        {
                            var connected = _simulation.FindConnectedParticles(i);
                            // Check if this group spans the X axis
                            bool spansX = CheckGroupSpansAxis(connected, 0);
                            if (spansX)
                            {
                                foreach (int idx in connected)
                                {
                                    if (!toRemove.Contains(idx))
                                        toRemove.Add(idx);
                                }
                            }
                        }
                    }

                    // Remove particles
                    foreach (int idx in toRemove)
                    {
                        _simulation.RemoveParticle(idx);
                    }

                    if (toRemove.Count > 0)
                    {
                        _score += ScorePerClear + (toRemove.Count * 10);
                        totalCleared += toRemove.Count;
                    }
                }

                // Also check Z axis
                if (_simulation.CheckSpansAxis(color, 2))
                {
                    var toRemove = new List<int>();

                    for (int i = 0; i < _simulation.Particles.Count; i++)
                    {
                        var p = _simulation.Particles[i];
                        if (p.IsActive && p.Color == color)
                        {
                            var connected = _simulation.FindConnectedParticles(i);
                            bool spansZ = CheckGroupSpansAxis(connected, 2);
                            if (spansZ)
                            {
                                foreach (int idx in connected)
                                {
                                    if (!toRemove.Contains(idx))
                                        toRemove.Add(idx);
                                }
                            }
                        }
                    }

                    foreach (int idx in toRemove)
                    {
                        _simulation.RemoveParticle(idx);
                    }

                    if (toRemove.Count > 0)
                    {
                        _score += ScorePerClear + (toRemove.Count * 10);
                        totalCleared += toRemove.Count;
                    }
                }
            }

            return totalCleared;
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

            double threshold = axis switch
            {
                0 => TerrariumWidth * 0.8,
                1 => TerrariumHeight * 0.8,
                2 => TerrariumDepth * 0.8,
                _ => 0
            };

            return (max - min) >= threshold;
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
        public void Reset()
        {
            _simulation.Clear();
            _score = 0;
            _turn = 0;
            _gameOver = false;
            _victory = false;
            _lastGroupId = 0;
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
