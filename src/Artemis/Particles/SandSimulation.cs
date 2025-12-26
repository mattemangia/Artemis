using System;
using System.Collections.Generic;
using Artemis.Core;

namespace Artemis.Particles
{
    /// <summary>
    /// Specialized particle simulation for sand and granular materials.
    /// Uses a grid-based approach for efficient particle-particle interactions.
    /// </summary>
    public class SandSimulation
    {
        #region Fields

        private readonly List<SandParticle> _particles;
        private readonly Dictionary<(int, int, int), List<int>> _grid;
        private readonly double _cellSize;

        #endregion

        #region Properties

        /// <summary>
        /// Gets all particles.
        /// </summary>
        public IReadOnlyList<SandParticle> Particles => _particles;

        /// <summary>
        /// Gets the number of particles.
        /// </summary>
        public int ParticleCount => _particles.Count;

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>
        /// Gets or sets the simulation bounds.
        /// </summary>
        public AABB Bounds { get; set; }

        /// <summary>
        /// Gets or sets the friction coefficient for particle sliding.
        /// </summary>
        public double Friction { get; set; } = 0.3;

        /// <summary>
        /// Gets or sets the restitution (bounciness) for collisions.
        /// </summary>
        public double Restitution { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets the default particle radius.
        /// </summary>
        public double ParticleRadius { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets the stacking angle in radians.
        /// Sand piles up at this angle (angle of repose).
        /// </summary>
        public double StackingAngle { get; set; } = Math.PI / 6; // 30 degrees

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new sand simulation.
        /// </summary>
        /// <param name="bounds">Simulation bounds.</param>
        /// <param name="cellSize">Grid cell size for spatial hashing.</param>
        public SandSimulation(AABB bounds, double cellSize = 0.5)
        {
            Bounds = bounds;
            _cellSize = cellSize;
            _particles = new List<SandParticle>();
            _grid = new Dictionary<(int, int, int), List<int>>();
        }

        #endregion

        #region Particle Management

        /// <summary>
        /// Adds a sand particle.
        /// </summary>
        /// <returns>The index of the added particle.</returns>
        public int AddParticle(Vector3D position, uint color, double? radius = null)
        {
            var particle = new SandParticle
            {
                Position = position,
                PreviousPosition = position,
                Velocity = Vector3D.Zero,
                Color = color,
                Radius = radius ?? ParticleRadius,
                IsActive = true,
                IsSettled = false,
                GroupId = -1
            };

            int index = _particles.Count;
            _particles.Add(particle);
            return index;
        }

        /// <summary>
        /// Adds a spherical cluster of sand particles.
        /// </summary>
        /// <returns>The starting index and count of added particles.</returns>
        public (int startIndex, int count) AddSandBall(
            Vector3D center,
            double ballRadius,
            uint color,
            double? particleRadius = null)
        {
            double pRadius = particleRadius ?? ParticleRadius;
            double spacing = pRadius * 2.1; // Slight gap to prevent initial overlap

            int startIndex = _particles.Count;
            int count = 0;

            // Fill sphere with particles in a grid pattern
            int steps = (int)Math.Ceiling(ballRadius / spacing);

            for (int x = -steps; x <= steps; x++)
            {
                for (int y = -steps; y <= steps; y++)
                {
                    for (int z = -steps; z <= steps; z++)
                    {
                        var offset = new Vector3D(x * spacing, y * spacing, z * spacing);
                        if (offset.Magnitude <= ballRadius - pRadius)
                        {
                            AddParticle(center + offset, color, pRadius);
                            count++;
                        }
                    }
                }
            }

            return (startIndex, count);
        }

        /// <summary>
        /// Removes a particle by index.
        /// </summary>
        public void RemoveParticle(int index)
        {
            if (index >= 0 && index < _particles.Count)
            {
                var p = _particles[index];
                p.IsActive = false;
                _particles[index] = p;
            }
        }

        /// <summary>
        /// Removes all particles with the specified group ID.
        /// </summary>
        public int RemoveGroup(int groupId)
        {
            int removed = 0;
            for (int i = 0; i < _particles.Count; i++)
            {
                if (_particles[i].GroupId == groupId && _particles[i].IsActive)
                {
                    var p = _particles[i];
                    p.IsActive = false;
                    _particles[i] = p;
                    removed++;
                }
            }
            return removed;
        }

        /// <summary>
        /// Clears all particles.
        /// </summary>
        public void Clear()
        {
            _particles.Clear();
            _grid.Clear();
        }

        #endregion

        #region Simulation

        /// <summary>
        /// Updates the simulation.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        /// <param name="subSteps">Number of sub-steps for stability.</param>
        public void Update(double deltaTime, int subSteps = 4)
        {
            double dt = deltaTime / subSteps;

            for (int step = 0; step < subSteps; step++)
            {
                // Rebuild spatial grid
                RebuildGrid();

                // Apply gravity and integrate
                for (int i = 0; i < _particles.Count; i++)
                {
                    if (!_particles[i].IsActive)
                        continue;

                    var p = _particles[i];

                    // Skip settled particles (optimization)
                    if (p.IsSettled)
                        continue;

                    // Apply gravity
                    p.Velocity += Gravity * dt;

                    // Verlet-style integration
                    p.PreviousPosition = p.Position;
                    p.Position += p.Velocity * dt;

                    _particles[i] = p;
                }

                // Handle collisions
                HandleParticleCollisions();
                HandleBoundCollisions();

                // Check for settling
                CheckSettling();
            }
        }

        private void RebuildGrid()
        {
            _grid.Clear();

            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive)
                    continue;

                var cell = GetCell(_particles[i].Position);
                if (!_grid.TryGetValue(cell, out var list))
                {
                    list = new List<int>();
                    _grid[cell] = list;
                }
                list.Add(i);
            }
        }

        private (int, int, int) GetCell(Vector3D position)
        {
            return (
                (int)Math.Floor(position.X / _cellSize),
                (int)Math.Floor(position.Y / _cellSize),
                (int)Math.Floor(position.Z / _cellSize)
            );
        }

        private void HandleParticleCollisions()
        {
            foreach (var (cell, indices) in _grid)
            {
                // Check within cell
                for (int i = 0; i < indices.Count; i++)
                {
                    for (int j = i + 1; j < indices.Count; j++)
                    {
                        ResolveCollision(indices[i], indices[j]);
                    }
                }

                // Check neighboring cells
                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            if (dx == 0 && dy == 0 && dz == 0)
                                continue;

                            var neighbor = (cell.Item1 + dx, cell.Item2 + dy, cell.Item3 + dz);
                            if (!_grid.TryGetValue(neighbor, out var neighborIndices))
                                continue;

                            foreach (var i in indices)
                            {
                                foreach (var j in neighborIndices)
                                {
                                    ResolveCollision(i, j);
                                }
                            }
                        }
                    }
                }
            }
        }

        private void ResolveCollision(int i, int j)
        {
            var a = _particles[i];
            var b = _particles[j];

            if (!a.IsActive || !b.IsActive)
                return;

            var delta = b.Position - a.Position;
            double distSq = delta.MagnitudeSquared;
            double minDist = a.Radius + b.Radius;

            if (distSq >= minDist * minDist || distSq < PhysicsConstants.Epsilon)
                return;

            double dist = Math.Sqrt(distSq);
            var normal = delta / dist;
            double penetration = minDist - dist;

            // Position correction
            var correction = normal * (penetration * 0.5);
            a.Position -= correction;
            b.Position += correction;

            // Velocity correction with friction
            var relVel = b.Velocity - a.Velocity;
            double velAlongNormal = Vector3D.Dot(relVel, normal);

            if (velAlongNormal > 0)
            {
                _particles[i] = a;
                _particles[j] = b;
                return;
            }

            double impulseMag = -(1 + Restitution) * velAlongNormal * 0.5;
            var impulse = normal * impulseMag;

            a.Velocity -= impulse;
            b.Velocity += impulse;

            // Tangential friction
            var tangent = relVel - normal * velAlongNormal;
            if (tangent.MagnitudeSquared > PhysicsConstants.Epsilon)
            {
                tangent = tangent.Normalized;
                double frictionImpulse = Friction * impulseMag;
                a.Velocity += tangent * frictionImpulse;
                b.Velocity -= tangent * frictionImpulse;
            }

            // Wake up settled particles
            a.IsSettled = false;
            b.IsSettled = false;

            _particles[i] = a;
            _particles[j] = b;
        }

        private void HandleBoundCollisions()
        {
            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive)
                    continue;

                var p = _particles[i];

                // Floor
                if (p.Position.Y - p.Radius < Bounds.Min.Y)
                {
                    p.Position.Y = Bounds.Min.Y + p.Radius;
                    p.Velocity.Y = -p.Velocity.Y * Restitution;
                    p.Velocity.X *= (1 - Friction);
                    p.Velocity.Z *= (1 - Friction);
                }

                // Ceiling
                if (p.Position.Y + p.Radius > Bounds.Max.Y)
                {
                    p.Position.Y = Bounds.Max.Y - p.Radius;
                    p.Velocity.Y = -p.Velocity.Y * Restitution;
                }

                // Walls
                if (p.Position.X - p.Radius < Bounds.Min.X)
                {
                    p.Position.X = Bounds.Min.X + p.Radius;
                    p.Velocity.X = -p.Velocity.X * Restitution;
                }
                if (p.Position.X + p.Radius > Bounds.Max.X)
                {
                    p.Position.X = Bounds.Max.X - p.Radius;
                    p.Velocity.X = -p.Velocity.X * Restitution;
                }
                if (p.Position.Z - p.Radius < Bounds.Min.Z)
                {
                    p.Position.Z = Bounds.Min.Z + p.Radius;
                    p.Velocity.Z = -p.Velocity.Z * Restitution;
                }
                if (p.Position.Z + p.Radius > Bounds.Max.Z)
                {
                    p.Position.Z = Bounds.Max.Z - p.Radius;
                    p.Velocity.Z = -p.Velocity.Z * Restitution;
                }

                _particles[i] = p;
            }
        }

        private void CheckSettling()
        {
            double settleThreshold = 0.01;
            double settleThresholdSq = settleThreshold * settleThreshold;

            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive || _particles[i].IsSettled)
                    continue;

                var p = _particles[i];
                if (p.Velocity.MagnitudeSquared < settleThresholdSq)
                {
                    p.Velocity = Vector3D.Zero;
                    p.IsSettled = true;
                    _particles[i] = p;
                }
            }
        }

        #endregion

        #region Queries

        /// <summary>
        /// Finds connected particles of the same color using flood fill.
        /// </summary>
        /// <param name="startIndex">Starting particle index.</param>
        /// <returns>List of connected particle indices.</returns>
        public List<int> FindConnectedParticles(int startIndex)
        {
            if (startIndex < 0 || startIndex >= _particles.Count || !_particles[startIndex].IsActive)
                return new List<int>();

            var result = new List<int>();
            var visited = new HashSet<int>();
            var queue = new Queue<int>();
            uint targetColor = _particles[startIndex].Color;

            queue.Enqueue(startIndex);
            visited.Add(startIndex);

            double connectionDistance = ParticleRadius * 2.5;
            double connectionDistSq = connectionDistance * connectionDistance;

            while (queue.Count > 0)
            {
                int current = queue.Dequeue();
                result.Add(current);

                var currentPos = _particles[current].Position;

                // Check nearby particles
                for (int i = 0; i < _particles.Count; i++)
                {
                    if (visited.Contains(i) || !_particles[i].IsActive)
                        continue;

                    if (_particles[i].Color != targetColor)
                        continue;

                    double distSq = Vector3D.DistanceSquared(currentPos, _particles[i].Position);
                    if (distSq <= connectionDistSq)
                    {
                        visited.Add(i);
                        queue.Enqueue(i);
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Checks if particles span from one side to another.
        /// </summary>
        /// <param name="color">The color to check.</param>
        /// <param name="axis">0=X, 1=Y, 2=Z</param>
        /// <returns>True if connected from min to max on the specified axis.</returns>
        public bool CheckSpansAxis(uint color, int axis)
        {
            // Find all particles touching the min side
            double minThreshold = GetAxisMin(axis) + ParticleRadius * 3;
            double maxThreshold = GetAxisMax(axis) - ParticleRadius * 3;

            var startParticles = new List<int>();
            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive || _particles[i].Color != color)
                    continue;

                double pos = GetAxisValue(_particles[i].Position, axis);
                if (pos <= minThreshold)
                    startParticles.Add(i);
            }

            // Flood fill from each starting particle
            foreach (int start in startParticles)
            {
                var connected = FindConnectedParticles(start);
                foreach (int idx in connected)
                {
                    double pos = GetAxisValue(_particles[idx].Position, axis);
                    if (pos >= maxThreshold)
                        return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Checks if there's a hole at the bottom of the terrarium.
        /// </summary>
        public bool HasHoleAtBottom()
        {
            double bottom = Bounds.Min.Y + ParticleRadius * 2;

            // Check if any position at the bottom is empty
            double checkRadius = ParticleRadius * 2;
            double stepX = ParticleRadius * 2;
            double stepZ = ParticleRadius * 2;

            for (double x = Bounds.Min.X + checkRadius; x < Bounds.Max.X - checkRadius; x += stepX)
            {
                for (double z = Bounds.Min.Z + checkRadius; z < Bounds.Max.Z - checkRadius; z += stepZ)
                {
                    var checkPos = new Vector3D(x, bottom, z);
                    bool hasParticle = false;

                    for (int i = 0; i < _particles.Count; i++)
                    {
                        if (!_particles[i].IsActive)
                            continue;

                        if (Vector3D.DistanceSquared(_particles[i].Position, checkPos) <
                            checkRadius * checkRadius)
                        {
                            hasParticle = true;
                            break;
                        }
                    }

                    if (!hasParticle)
                        return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Checks if the terrarium is full (particles reaching the top).
        /// </summary>
        public bool IsFull()
        {
            double topThreshold = Bounds.Max.Y - ParticleRadius * 4;

            for (int i = 0; i < _particles.Count; i++)
            {
                if (_particles[i].IsActive && _particles[i].Position.Y >= topThreshold)
                    return true;
            }

            return false;
        }

        private double GetAxisValue(Vector3D v, int axis) => axis switch
        {
            0 => v.X,
            1 => v.Y,
            2 => v.Z,
            _ => 0
        };

        private double GetAxisMin(int axis) => axis switch
        {
            0 => Bounds.Min.X,
            1 => Bounds.Min.Y,
            2 => Bounds.Min.Z,
            _ => 0
        };

        private double GetAxisMax(int axis) => axis switch
        {
            0 => Bounds.Max.X,
            1 => Bounds.Max.Y,
            2 => Bounds.Max.Z,
            _ => 0
        };

        #endregion
    }

    /// <summary>
    /// Represents a sand particle.
    /// </summary>
    public struct SandParticle
    {
        public Vector3D Position;
        public Vector3D PreviousPosition;
        public Vector3D Velocity;
        public double Radius;
        public uint Color;
        public bool IsActive;
        public bool IsSettled;
        public int GroupId;
    }
}
