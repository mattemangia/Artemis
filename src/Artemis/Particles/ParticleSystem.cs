using System;
using System.Collections.Generic;
using Artemis.Core;
using Artemis.Forces;
using Artemis.Materials;

namespace Artemis.Particles
{
    /// <summary>
    /// Manages a collection of particles with efficient batch updates.
    /// </summary>
    public class ParticleSystem
    {
        #region Fields

        private Particle[] _particles;
        private int _aliveCount;
        private readonly List<IForce> _forces;
        private readonly Random _random;

        #endregion

        #region Properties

        /// <summary>
        /// Gets or sets the maximum number of particles.
        /// </summary>
        public int MaxParticles { get; private set; }

        /// <summary>
        /// Gets the number of alive particles.
        /// </summary>
        public int AliveCount => _aliveCount;

        /// <summary>
        /// Gets the particle array (for direct access in renderers).
        /// </summary>
        public Particle[] Particles => _particles;

        /// <summary>
        /// Gets or sets the default particle lifetime.
        /// </summary>
        public double DefaultLifetime { get; set; } = 5.0;

        /// <summary>
        /// Gets or sets the default particle radius.
        /// </summary>
        public double DefaultRadius { get; set; } = 0.1;

        /// <summary>
        /// Gets or sets the default particle mass.
        /// </summary>
        public double DefaultMass { get; set; } = 1.0;

        /// <summary>
        /// Gets or sets the damping coefficient (0-1).
        /// </summary>
        public double Damping { get; set; } = 0.99;

        /// <summary>
        /// Gets or sets whether to use Verlet integration.
        /// </summary>
        public bool UseVerletIntegration { get; set; } = false;

        /// <summary>
        /// Gets or sets the collision restitution for particle-particle collisions.
        /// </summary>
        public double CollisionRestitution { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>
        /// Gets or sets the bounding box for the particle system.
        /// Particles outside this box can be killed or bounced.
        /// </summary>
        public AABB? Bounds { get; set; }

        /// <summary>
        /// Gets or sets whether particles bounce off bounds.
        /// </summary>
        public bool BounceOffBounds { get; set; } = true;

        /// <summary>
        /// Gets or sets the bound bounce restitution.
        /// </summary>
        public double BoundRestitution { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets whether particle-particle collisions are enabled.
        /// </summary>
        public bool EnableParticleCollisions { get; set; } = false;

        /// <summary>
        /// Event raised when a particle dies.
        /// </summary>
        public event Action<Particle>? OnParticleDeath;

        /// <summary>
        /// Event raised when a particle collides with bounds.
        /// </summary>
        public event Action<int, Vector3D>? OnBoundCollision;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new particle system with the specified capacity.
        /// </summary>
        /// <param name="maxParticles">Maximum number of particles.</param>
        public ParticleSystem(int maxParticles = 10000)
        {
            MaxParticles = maxParticles;
            _particles = new Particle[maxParticles];
            _forces = new List<IForce>();
            _random = new Random();
            _aliveCount = 0;
        }

        #endregion

        #region Force Management

        /// <summary>
        /// Adds a force to the particle system.
        /// </summary>
        public void AddForce(IForce force)
        {
            _forces.Add(force);
        }

        /// <summary>
        /// Removes a force from the particle system.
        /// </summary>
        public bool RemoveForce(IForce force)
        {
            return _forces.Remove(force);
        }

        /// <summary>
        /// Clears all forces.
        /// </summary>
        public void ClearForces()
        {
            _forces.Clear();
        }

        #endregion

        #region Particle Spawning

        /// <summary>
        /// Spawns a single particle.
        /// </summary>
        /// <returns>The index of the spawned particle, or -1 if system is full.</returns>
        public int Spawn(
            Vector3D position,
            Vector3D velocity,
            double? mass = null,
            double? radius = null,
            double? lifetime = null,
            uint color = 0xFFFFFFFF,
            ParticleFlags flags = ParticleFlags.AffectedByGravity | ParticleFlags.CollidesWithWorld)
        {
            // Find a dead particle slot
            int index = FindDeadSlot();
            if (index < 0)
                return -1;

            _particles[index] = Particle.Create(
                position,
                velocity,
                mass ?? DefaultMass,
                radius ?? DefaultRadius,
                lifetime ?? DefaultLifetime,
                color
            );
            _particles[index].Flags = flags;
            _aliveCount++;

            return index;
        }

        /// <summary>
        /// Spawns multiple particles at once.
        /// </summary>
        public int SpawnBurst(
            Vector3D position,
            int count,
            double spreadAngle = Math.PI * 2,
            double minSpeed = 1.0,
            double maxSpeed = 5.0,
            double? mass = null,
            double? radius = null,
            double? lifetime = null,
            uint color = 0xFFFFFFFF)
        {
            int spawned = 0;

            for (int i = 0; i < count; i++)
            {
                // Random direction on sphere
                double theta = _random.NextDouble() * spreadAngle;
                double phi = Math.Acos(2 * _random.NextDouble() - 1);
                double speed = minSpeed + _random.NextDouble() * (maxSpeed - minSpeed);

                var velocity = new Vector3D(
                    Math.Sin(phi) * Math.Cos(theta) * speed,
                    Math.Sin(phi) * Math.Sin(theta) * speed,
                    Math.Cos(phi) * speed
                );

                if (Spawn(position, velocity, mass, radius, lifetime, color) >= 0)
                    spawned++;
            }

            return spawned;
        }

        /// <summary>
        /// Spawns particles in a sphere shape.
        /// </summary>
        public int SpawnSphere(
            Vector3D center,
            double sphereRadius,
            int count,
            double? particleMass = null,
            double? particleRadius = null,
            uint color = 0xFFFFFFFF)
        {
            int spawned = 0;

            for (int i = 0; i < count; i++)
            {
                // Random point in sphere using rejection sampling
                Vector3D offset;
                do
                {
                    offset = new Vector3D(
                        (_random.NextDouble() * 2 - 1) * sphereRadius,
                        (_random.NextDouble() * 2 - 1) * sphereRadius,
                        (_random.NextDouble() * 2 - 1) * sphereRadius
                    );
                } while (offset.MagnitudeSquared > sphereRadius * sphereRadius);

                if (Spawn(center + offset, Vector3D.Zero, particleMass, particleRadius, null, color) >= 0)
                    spawned++;
            }

            return spawned;
        }

        private int FindDeadSlot()
        {
            for (int i = 0; i < MaxParticles; i++)
            {
                if (!_particles[i].IsAlive)
                    return i;
            }
            return -1;
        }

        #endregion

        #region Update

        /// <summary>
        /// Updates all particles.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void Update(double deltaTime)
        {
            _aliveCount = 0;

            // Apply forces and integrate
            for (int i = 0; i < MaxParticles; i++)
            {
                if (!_particles[i].IsAlive)
                    continue;

                // Apply gravity if flagged
                if (_particles[i].Flags.HasFlag(ParticleFlags.AffectedByGravity))
                {
                    _particles[i].ApplyForce(Gravity * _particles[i].Mass);
                }

                // Apply external forces
                foreach (var force in _forces)
                {
                    if (force.Enabled)
                    {
                        var f = force.Calculate(
                            _particles[i].Position,
                            _particles[i].Velocity,
                            _particles[i].Mass
                        );
                        _particles[i].ApplyForce(f);
                    }
                }

                // Integrate
                if (UseVerletIntegration)
                    _particles[i].IntegrateVerlet(deltaTime, Damping);
                else
                    _particles[i].IntegrateEuler(deltaTime, Damping);

                // Handle bounds
                if (Bounds.HasValue && _particles[i].Flags.HasFlag(ParticleFlags.CollidesWithWorld))
                {
                    HandleBoundsCollision(ref _particles[i], i);
                }

                // Count alive
                if (_particles[i].IsAlive)
                    _aliveCount++;
                else
                    OnParticleDeath?.Invoke(_particles[i]);
            }

            // Handle particle-particle collisions
            if (EnableParticleCollisions)
            {
                HandleParticleCollisions();
            }
        }

        private void HandleBoundsCollision(ref Particle particle, int index)
        {
            var bounds = Bounds!.Value;
            bool collided = false;
            var normal = Vector3D.Zero;

            // Check each axis
            if (particle.Position.X - particle.Radius < bounds.Min.X)
            {
                particle.Position.X = bounds.Min.X + particle.Radius;
                particle.Velocity.X = -particle.Velocity.X * BoundRestitution;
                normal = Vector3D.Right;
                collided = true;
            }
            else if (particle.Position.X + particle.Radius > bounds.Max.X)
            {
                particle.Position.X = bounds.Max.X - particle.Radius;
                particle.Velocity.X = -particle.Velocity.X * BoundRestitution;
                normal = Vector3D.Left;
                collided = true;
            }

            if (particle.Position.Y - particle.Radius < bounds.Min.Y)
            {
                particle.Position.Y = bounds.Min.Y + particle.Radius;
                particle.Velocity.Y = -particle.Velocity.Y * BoundRestitution;
                normal = Vector3D.Up;
                collided = true;
            }
            else if (particle.Position.Y + particle.Radius > bounds.Max.Y)
            {
                particle.Position.Y = bounds.Max.Y - particle.Radius;
                particle.Velocity.Y = -particle.Velocity.Y * BoundRestitution;
                normal = Vector3D.Down;
                collided = true;
            }

            if (particle.Position.Z - particle.Radius < bounds.Min.Z)
            {
                particle.Position.Z = bounds.Min.Z + particle.Radius;
                particle.Velocity.Z = -particle.Velocity.Z * BoundRestitution;
                normal = Vector3D.Forward;
                collided = true;
            }
            else if (particle.Position.Z + particle.Radius > bounds.Max.Z)
            {
                particle.Position.Z = bounds.Max.Z - particle.Radius;
                particle.Velocity.Z = -particle.Velocity.Z * BoundRestitution;
                normal = Vector3D.Backward;
                collided = true;
            }

            if (collided)
            {
                OnBoundCollision?.Invoke(index, normal);

                if (!BounceOffBounds)
                    particle.Kill();
            }
        }

        private void HandleParticleCollisions()
        {
            // Simple O(n²) collision detection - for large systems, use spatial hashing
            for (int i = 0; i < MaxParticles; i++)
            {
                if (!_particles[i].IsAlive || !_particles[i].Flags.HasFlag(ParticleFlags.CollidesWithParticles))
                    continue;

                for (int j = i + 1; j < MaxParticles; j++)
                {
                    if (!_particles[j].IsAlive || !_particles[j].Flags.HasFlag(ParticleFlags.CollidesWithParticles))
                        continue;

                    ResolveParticleCollision(ref _particles[i], ref _particles[j]);
                }
            }
        }

        private void ResolveParticleCollision(ref Particle a, ref Particle b)
        {
            var delta = b.Position - a.Position;
            double distSq = delta.MagnitudeSquared;
            double radiusSum = a.Radius + b.Radius;

            if (distSq >= radiusSum * radiusSum)
                return;

            double dist = Math.Sqrt(distSq);
            if (dist < PhysicsConstants.Epsilon)
            {
                // Particles at same position, separate randomly
                delta = new Vector3D(
                    _random.NextDouble() - 0.5,
                    _random.NextDouble() - 0.5,
                    _random.NextDouble() - 0.5
                ).Normalized;
                dist = 0.001;
            }

            var normal = delta / dist;
            double penetration = radiusSum - dist;

            // Separate particles
            double totalMass = a.InverseMass + b.InverseMass;
            if (totalMass > 0)
            {
                var correction = normal * (penetration / totalMass);
                a.Position -= correction * a.InverseMass;
                b.Position += correction * b.InverseMass;
            }

            // Impulse-based velocity resolution
            var relVel = b.Velocity - a.Velocity;
            double velAlongNormal = Vector3D.Dot(relVel, normal);

            if (velAlongNormal > 0)
                return; // Separating

            double j = -(1 + CollisionRestitution) * velAlongNormal / totalMass;
            var impulse = normal * j;

            a.Velocity -= impulse * a.InverseMass;
            b.Velocity += impulse * b.InverseMass;
        }

        #endregion

        #region Utility Methods

        /// <summary>
        /// Gets a particle by index.
        /// </summary>
        public ref Particle GetParticle(int index) => ref _particles[index];

        /// <summary>
        /// Clears all particles.
        /// </summary>
        public void Clear()
        {
            for (int i = 0; i < MaxParticles; i++)
            {
                _particles[i].IsAlive = false;
            }
            _aliveCount = 0;
        }

        /// <summary>
        /// Kills all particles.
        /// </summary>
        public void KillAll()
        {
            Clear();
        }

        /// <summary>
        /// Applies an explosion force to all particles.
        /// </summary>
        public void ApplyExplosion(Vector3D center, double force, double radius)
        {
            double radiusSq = radius * radius;

            for (int i = 0; i < MaxParticles; i++)
            {
                if (!_particles[i].IsAlive)
                    continue;

                var delta = _particles[i].Position - center;
                double distSq = delta.MagnitudeSquared;

                if (distSq < radiusSq && distSq > PhysicsConstants.Epsilon)
                {
                    double dist = Math.Sqrt(distSq);
                    double falloff = 1.0 - dist / radius;
                    var direction = delta / dist;
                    _particles[i].ApplyForce(direction * force * falloff);
                }
            }
        }

        /// <summary>
        /// Gets all particles of a specific color.
        /// </summary>
        public List<int> GetParticlesByColor(uint color)
        {
            var result = new List<int>();
            for (int i = 0; i < MaxParticles; i++)
            {
                if (_particles[i].IsAlive && _particles[i].Color == color)
                    result.Add(i);
            }
            return result;
        }

        /// <summary>
        /// Resizes the particle system.
        /// </summary>
        public void Resize(int newMaxParticles)
        {
            if (newMaxParticles == MaxParticles)
                return;

            var newParticles = new Particle[newMaxParticles];
            int copyCount = Math.Min(MaxParticles, newMaxParticles);

            Array.Copy(_particles, newParticles, copyCount);

            _particles = newParticles;
            MaxParticles = newMaxParticles;
        }

        #endregion
    }
}
