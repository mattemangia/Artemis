using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Collision;
using Artemis.Core;
using Artemis.Forces;
using Artemis.Particles;

namespace Artemis.Simulation
{
    /// <summary>
    /// The main physics simulation world that manages bodies, forces, and updates.
    /// </summary>
    public class PhysicsWorld
    {
        #region Fields

        private readonly List<IPhysicsBody> _bodies;
        private readonly List<IForce> _globalForces;
        private readonly List<ParticleSystem> _particleSystems;
        private readonly List<CollisionInfo> _collisions;
        private double _accumulator;

        #endregion

        #region Properties

        /// <summary>
        /// Gets all bodies in the world.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> Bodies => _bodies;

        /// <summary>
        /// Gets all global forces.
        /// </summary>
        public IReadOnlyList<IForce> GlobalForces => _globalForces;

        /// <summary>
        /// Gets all particle systems.
        /// </summary>
        public IReadOnlyList<ParticleSystem> ParticleSystems => _particleSystems;

        /// <summary>
        /// Gets the collisions from the last update.
        /// </summary>
        public IReadOnlyList<CollisionInfo> LastCollisions => _collisions;

        /// <summary>
        /// Gets or sets the fixed time step for physics simulation.
        /// </summary>
        public double FixedTimeStep { get; set; } = PhysicsConstants.DefaultFixedTimeStep;

        /// <summary>
        /// Gets or sets the maximum number of sub-steps per frame.
        /// </summary>
        public int MaxSubSteps { get; set; } = PhysicsConstants.DefaultMaxSubSteps;

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity { get; set; } = new(0, -PhysicsConstants.EarthGravity, 0);

        /// <summary>
        /// Gets or sets whether physics simulation is paused.
        /// </summary>
        public bool IsPaused { get; set; }

        /// <summary>
        /// Gets or sets the velocity threshold for collision restitution.
        /// </summary>
        public double RestitutionThreshold { get; set; } = 1.0;

        /// <summary>
        /// Gets or sets the world bounds. Bodies outside may be deactivated.
        /// </summary>
        public AABB? WorldBounds { get; set; }

        /// <summary>
        /// Gets or sets whether to use continuous collision detection.
        /// </summary>
        public bool UseCCD { get; set; } = false;

        /// <summary>
        /// Gets the total simulation time.
        /// </summary>
        public double TotalTime { get; private set; }

        /// <summary>
        /// Gets the number of physics steps taken.
        /// </summary>
        public long StepCount { get; private set; }

        /// <summary>
        /// Event raised when a collision occurs.
        /// </summary>
        public event Action<CollisionInfo>? OnCollision;

        /// <summary>
        /// Event raised before each physics step.
        /// </summary>
        public event Action<double>? OnPreStep;

        /// <summary>
        /// Event raised after each physics step.
        /// </summary>
        public event Action<double>? OnPostStep;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new physics world with default settings.
        /// </summary>
        public PhysicsWorld()
        {
            _bodies = new List<IPhysicsBody>();
            _globalForces = new List<IForce>();
            _particleSystems = new List<ParticleSystem>();
            _collisions = new List<CollisionInfo>();
        }

        /// <summary>
        /// Creates a new physics world with custom gravity.
        /// </summary>
        /// <param name="gravity">The gravity vector.</param>
        public PhysicsWorld(Vector3D gravity) : this()
        {
            Gravity = gravity;
        }

        #endregion

        #region Body Management

        /// <summary>
        /// Adds a body to the world.
        /// </summary>
        public void AddBody(IPhysicsBody body)
        {
            if (!_bodies.Contains(body))
            {
                _bodies.Add(body);
            }
        }

        /// <summary>
        /// Removes a body from the world.
        /// </summary>
        public bool RemoveBody(IPhysicsBody body)
        {
            return _bodies.Remove(body);
        }

        /// <summary>
        /// Gets a body by its ID.
        /// </summary>
        public IPhysicsBody? GetBody(string id)
        {
            return _bodies.Find(b => b.Id == id);
        }

        /// <summary>
        /// Clears all bodies from the world.
        /// </summary>
        public void ClearBodies()
        {
            _bodies.Clear();
        }

        #endregion

        #region Force Management

        /// <summary>
        /// Adds a global force that affects all bodies.
        /// </summary>
        public void AddForce(IForce force)
        {
            if (!_globalForces.Contains(force))
            {
                _globalForces.Add(force);
            }
        }

        /// <summary>
        /// Removes a global force.
        /// </summary>
        public bool RemoveForce(IForce force)
        {
            return _globalForces.Remove(force);
        }

        /// <summary>
        /// Clears all global forces.
        /// </summary>
        public void ClearForces()
        {
            _globalForces.Clear();
        }

        #endregion

        #region Particle System Management

        /// <summary>
        /// Adds a particle system to the world.
        /// </summary>
        public void AddParticleSystem(ParticleSystem system)
        {
            if (!_particleSystems.Contains(system))
            {
                _particleSystems.Add(system);
            }
        }

        /// <summary>
        /// Removes a particle system from the world.
        /// </summary>
        public bool RemoveParticleSystem(ParticleSystem system)
        {
            return _particleSystems.Remove(system);
        }

        #endregion

        #region Simulation

        /// <summary>
        /// Updates the physics simulation.
        /// </summary>
        /// <param name="deltaTime">Time since last update in seconds.</param>
        public void Update(double deltaTime)
        {
            if (IsPaused)
                return;

            // Accumulate time
            _accumulator += deltaTime;

            // Limit accumulator to prevent spiral of death
            double maxAccumulator = FixedTimeStep * MaxSubSteps;
            if (_accumulator > maxAccumulator)
                _accumulator = maxAccumulator;

            // Fixed time step simulation
            while (_accumulator >= FixedTimeStep)
            {
                Step(FixedTimeStep);
                _accumulator -= FixedTimeStep;
            }

            // Update particle systems with remaining time
            foreach (var system in _particleSystems)
            {
                system.Update(deltaTime);
            }
        }

        /// <summary>
        /// Performs a single physics step.
        /// </summary>
        /// <param name="dt">Time step in seconds.</param>
        public void Step(double dt)
        {
            OnPreStep?.Invoke(dt);

            // Clear previous collision info
            _collisions.Clear();

            // Apply gravity and global forces to all bodies
            ApplyForces(dt);

            // Integrate velocities and positions
            IntegrateBodies(dt);

            // Detect collisions
            DetectCollisions();

            // Resolve collisions
            ResolveCollisions();

            // Update simulation time
            TotalTime += dt;
            StepCount++;

            OnPostStep?.Invoke(dt);
        }

        private void ApplyForces(double dt)
        {
            foreach (var body in _bodies)
            {
                if (body.BodyType != BodyType.Dynamic || !body.IsActive || body.IsSleeping)
                    continue;

                // Apply gravity
                body.ApplyForce(Gravity * body.Mass);

                // Apply global forces
                foreach (var force in _globalForces)
                {
                    if (force.Enabled)
                    {
                        var f = force.Calculate(body.Position, body.Velocity, body.Mass);
                        body.ApplyForce(f);
                    }
                }
            }
        }

        private void IntegrateBodies(double dt)
        {
            foreach (var body in _bodies)
            {
                if (!body.IsActive)
                    continue;

                body.Integrate(dt);

                // Check world bounds
                if (WorldBounds.HasValue && !WorldBounds.Value.Contains(body.Position))
                {
                    // Option: deactivate or wrap around
                    // body.IsActive = false;
                }
            }
        }

        private void DetectCollisions()
        {
            // Broad phase: simple O(n²) for now
            // For larger simulations, implement spatial partitioning (octree, grid)
            for (int i = 0; i < _bodies.Count; i++)
            {
                var bodyA = _bodies[i];
                if (!bodyA.IsActive)
                    continue;

                for (int j = i + 1; j < _bodies.Count; j++)
                {
                    var bodyB = _bodies[j];
                    if (!bodyB.IsActive)
                        continue;

                    // Skip if both are static
                    if (bodyA.BodyType == BodyType.Static && bodyB.BodyType == BodyType.Static)
                        continue;

                    var collision = CollisionDetector.Detect(bodyA, bodyB);
                    if (collision.HasCollision)
                    {
                        _collisions.Add(collision);
                    }
                }
            }
        }

        private void ResolveCollisions()
        {
            foreach (var collision in _collisions)
            {
                CollisionResolver.Resolve(collision, RestitutionThreshold);
                OnCollision?.Invoke(collision);
            }
        }

        #endregion

        #region Queries

        /// <summary>
        /// Casts a ray through the world.
        /// </summary>
        /// <param name="origin">Ray origin.</param>
        /// <param name="direction">Ray direction (normalized).</param>
        /// <param name="maxDistance">Maximum ray distance.</param>
        /// <returns>Hit info or null if no hit.</returns>
        public RaycastHit? Raycast(Vector3D origin, Vector3D direction, double maxDistance = 1000)
        {
            RaycastHit? closest = null;
            double closestDist = maxDistance;

            foreach (var body in _bodies)
            {
                if (!body.IsActive)
                    continue;

                if (body is RigidBody rb)
                {
                    if (CollisionDetector.Raycast(origin, direction, rb, closestDist,
                        out double dist, out Vector3D point, out Vector3D normal))
                    {
                        if (dist < closestDist)
                        {
                            closestDist = dist;
                            closest = new RaycastHit
                            {
                                Body = body,
                                Distance = dist,
                                Point = point,
                                Normal = normal
                            };
                        }
                    }
                }
            }

            return closest;
        }

        /// <summary>
        /// Gets all bodies within a sphere.
        /// </summary>
        public List<IPhysicsBody> QuerySphere(Vector3D center, double radius)
        {
            var results = new List<IPhysicsBody>();
            double radiusSq = radius * radius;

            foreach (var body in _bodies)
            {
                if (!body.IsActive)
                    continue;

                if (Vector3D.DistanceSquared(center, body.Position) <= radiusSq)
                {
                    results.Add(body);
                }
            }

            return results;
        }

        /// <summary>
        /// Gets all bodies within an AABB.
        /// </summary>
        public List<IPhysicsBody> QueryAABB(AABB bounds)
        {
            var results = new List<IPhysicsBody>();

            foreach (var body in _bodies)
            {
                if (!body.IsActive)
                    continue;

                if (bounds.Intersects(body.BoundingBox))
                {
                    results.Add(body);
                }
            }

            return results;
        }

        #endregion

        #region Utility

        /// <summary>
        /// Resets the simulation.
        /// </summary>
        public void Reset()
        {
            TotalTime = 0;
            StepCount = 0;
            _accumulator = 0;
            _collisions.Clear();

            // Wake up all bodies
            foreach (var body in _bodies)
            {
                body.WakeUp();
            }
        }

        /// <summary>
        /// Creates a static ground plane.
        /// </summary>
        public RigidBody CreateGround(double width = 100, double depth = 100, double height = 1)
        {
            var ground = RigidBody.CreateStaticBox(
                new Vector3D(0, -height / 2, 0),
                new Vector3D(width / 2, height / 2, depth / 2)
            );
            AddBody(ground);
            return ground;
        }

        /// <summary>
        /// Creates walls around an area.
        /// </summary>
        public void CreateWalls(double width, double height, double depth, double thickness = 1)
        {
            // Left wall
            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(-width / 2 - thickness / 2, height / 2, 0),
                new Vector3D(thickness / 2, height / 2, depth / 2)
            ));

            // Right wall
            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(width / 2 + thickness / 2, height / 2, 0),
                new Vector3D(thickness / 2, height / 2, depth / 2)
            ));

            // Front wall
            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(0, height / 2, depth / 2 + thickness / 2),
                new Vector3D(width / 2, height / 2, thickness / 2)
            ));

            // Back wall
            AddBody(RigidBody.CreateStaticBox(
                new Vector3D(0, height / 2, -depth / 2 - thickness / 2),
                new Vector3D(width / 2, height / 2, thickness / 2)
            ));
        }

        #endregion
    }

    /// <summary>
    /// Information about a raycast hit.
    /// </summary>
    public struct RaycastHit
    {
        /// <summary>The body that was hit.</summary>
        public IPhysicsBody Body;

        /// <summary>Distance from ray origin to hit point.</summary>
        public double Distance;

        /// <summary>Hit point in world coordinates.</summary>
        public Vector3D Point;

        /// <summary>Surface normal at hit point.</summary>
        public Vector3D Normal;
    }
}
