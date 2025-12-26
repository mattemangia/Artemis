using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;
using Artemis.Forces;

namespace Artemis.Particles
{
    /// <summary>
    /// Smoothed Particle Hydrodynamics (SPH) fluid simulation.
    /// Simulates realistic fluid behavior with pressure, viscosity, and surface tension.
    /// </summary>
    public class FluidSimulation
    {
        #region Fields

        private readonly List<FluidParticle> _particles;
        private readonly Dictionary<(int, int, int), List<int>> _grid;
        private readonly List<IPhysicsBody> _colliders;
        private readonly List<IForce> _externalForces;
        private readonly double _cellSize;

        #endregion

        #region Properties

        /// <summary>
        /// Gets all particles.
        /// </summary>
        public IReadOnlyList<FluidParticle> Particles => _particles;

        /// <summary>
        /// Gets the number of particles.
        /// </summary>
        public int ParticleCount => _particles.Count;

        /// <summary>
        /// Gets the colliders list.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> Colliders => _colliders;

        /// <summary>
        /// Gets or sets the simulation bounds.
        /// </summary>
        public AABB Bounds { get; set; }

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>
        /// Gets or sets the rest density of the fluid in kg/m³.
        /// Water = 1000, Oil = 900, Honey = 1400
        /// </summary>
        public double RestDensity { get; set; } = 1000.0;

        /// <summary>
        /// Gets or sets the gas constant for pressure calculation.
        /// Higher = stiffer fluid.
        /// </summary>
        public double GasConstant { get; set; } = 2000.0;

        /// <summary>
        /// Gets or sets the viscosity coefficient.
        /// Water = 0.001, Oil = 0.1, Honey = 10
        /// </summary>
        public double Viscosity { get; set; } = 0.001;

        /// <summary>
        /// Gets or sets the surface tension coefficient.
        /// </summary>
        public double SurfaceTension { get; set; } = 0.0728;

        /// <summary>
        /// Gets or sets the smoothing radius (interaction distance).
        /// </summary>
        public double SmoothingRadius { get; set; } = 0.2;

        /// <summary>
        /// Gets or sets the particle mass.
        /// </summary>
        public double ParticleMass { get; set; } = 0.02;

        /// <summary>
        /// Gets or sets the boundary restitution.
        /// </summary>
        public double BoundaryRestitution { get; set; } = 0.3;

        /// <summary>
        /// Gets or sets the collision restitution with rigid bodies.
        /// </summary>
        public double CollisionRestitution { get; set; } = 0.2;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new fluid simulation.
        /// </summary>
        public FluidSimulation(AABB bounds, double smoothingRadius = 0.2)
        {
            Bounds = bounds;
            SmoothingRadius = smoothingRadius;
            _cellSize = smoothingRadius;
            _particles = new List<FluidParticle>();
            _grid = new Dictionary<(int, int, int), List<int>>();
            _colliders = new List<IPhysicsBody>();
            _externalForces = new List<IForce>();
        }

        #endregion

        #region Particle Management

        /// <summary>
        /// Adds a fluid particle.
        /// </summary>
        public int AddParticle(Vector3D position, Vector3D velocity = default)
        {
            var particle = new FluidParticle
            {
                Position = position,
                Velocity = velocity,
                Force = Vector3D.Zero,
                Density = RestDensity,
                Pressure = 0,
                IsActive = true
            };

            int index = _particles.Count;
            _particles.Add(particle);
            return index;
        }

        /// <summary>
        /// Adds a block of fluid particles.
        /// </summary>
        public int AddFluidBlock(Vector3D min, Vector3D max, double spacing)
        {
            int count = 0;
            double halfSpacing = spacing * 0.5;

            for (double x = min.X + halfSpacing; x < max.X; x += spacing)
            {
                for (double y = min.Y + halfSpacing; y < max.Y; y += spacing)
                {
                    for (double z = min.Z + halfSpacing; z < max.Z; z += spacing)
                    {
                        AddParticle(new Vector3D(x, y, z));
                        count++;
                    }
                }
            }

            return count;
        }

        /// <summary>
        /// Adds a sphere of fluid particles.
        /// </summary>
        public int AddFluidSphere(Vector3D center, double radius, double spacing)
        {
            int count = 0;
            double radiusSq = radius * radius;

            for (double x = center.X - radius; x <= center.X + radius; x += spacing)
            {
                for (double y = center.Y - radius; y <= center.Y + radius; y += spacing)
                {
                    for (double z = center.Z - radius; z <= center.Z + radius; z += spacing)
                    {
                        var pos = new Vector3D(x, y, z);
                        if (Vector3D.DistanceSquared(pos, center) <= radiusSq)
                        {
                            AddParticle(pos);
                            count++;
                        }
                    }
                }
            }

            return count;
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

        #region Collider Management

        /// <summary>
        /// Adds a rigid body as a collider.
        /// </summary>
        public void AddCollider(IPhysicsBody body)
        {
            if (!_colliders.Contains(body))
                _colliders.Add(body);
        }

        /// <summary>
        /// Removes a collider.
        /// </summary>
        public bool RemoveCollider(IPhysicsBody body)
        {
            return _colliders.Remove(body);
        }

        /// <summary>
        /// Clears all colliders.
        /// </summary>
        public void ClearColliders()
        {
            _colliders.Clear();
        }

        #endregion

        #region Force Management

        /// <summary>
        /// Adds an external force affecting the fluid.
        /// </summary>
        public void AddForce(IForce force)
        {
            if (!_externalForces.Contains(force))
                _externalForces.Add(force);
        }

        /// <summary>
        /// Removes an external force.
        /// </summary>
        public bool RemoveForce(IForce force)
        {
            return _externalForces.Remove(force);
        }

        #endregion

        #region Simulation

        /// <summary>
        /// Updates the fluid simulation.
        /// </summary>
        public void Update(double deltaTime, int subSteps = 2)
        {
            double dt = deltaTime / subSteps;

            for (int step = 0; step < subSteps; step++)
            {
                // Build spatial grid
                RebuildGrid();

                // Compute densities and pressures
                ComputeDensityPressure();

                // Compute forces
                ComputeForces();

                // Apply external forces
                ApplyExternalForces();

                // Integrate
                Integrate(dt);

                // Handle collisions
                HandleBoundaryCollisions();
                HandleColliderCollisions();
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

        private void ComputeDensityPressure()
        {
            double h = SmoothingRadius;
            double h2 = h * h;
            double poly6Coeff = 315.0 / (64.0 * Math.PI * Math.Pow(h, 9));

            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive)
                    continue;

                var pi = _particles[i];
                double density = 0;

                // Get neighboring cells
                var cell = GetCell(pi.Position);
                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            var neighborCell = (cell.Item1 + dx, cell.Item2 + dy, cell.Item3 + dz);
                            if (!_grid.TryGetValue(neighborCell, out var neighbors))
                                continue;

                            foreach (int j in neighbors)
                            {
                                var pj = _particles[j];
                                if (!pj.IsActive)
                                    continue;

                                double r2 = Vector3D.DistanceSquared(pi.Position, pj.Position);
                                if (r2 < h2)
                                {
                                    // Poly6 kernel
                                    double diff = h2 - r2;
                                    density += ParticleMass * poly6Coeff * diff * diff * diff;
                                }
                            }
                        }
                    }
                }

                pi.Density = density;
                // Tait equation for pressure
                pi.Pressure = GasConstant * (density - RestDensity);

                _particles[i] = pi;
            }
        }

        private void ComputeForces()
        {
            double h = SmoothingRadius;
            double h2 = h * h;
            double spikyCoeff = -45.0 / (Math.PI * Math.Pow(h, 6));
            double viscosityCoeff = 45.0 / (Math.PI * Math.Pow(h, 6));

            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive)
                    continue;

                var pi = _particles[i];
                var pressureForce = Vector3D.Zero;
                var viscosityForce = Vector3D.Zero;

                var cell = GetCell(pi.Position);
                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            var neighborCell = (cell.Item1 + dx, cell.Item2 + dy, cell.Item3 + dz);
                            if (!_grid.TryGetValue(neighborCell, out var neighbors))
                                continue;

                            foreach (int j in neighbors)
                            {
                                if (i == j)
                                    continue;

                                var pj = _particles[j];
                                if (!pj.IsActive)
                                    continue;

                                var rij = pi.Position - pj.Position;
                                double r = rij.Magnitude;

                                if (r < h && r > PhysicsConstants.Epsilon)
                                {
                                    var rijNorm = rij / r;

                                    // Pressure force (Spiky gradient)
                                    double pressureGrad = spikyCoeff * (h - r) * (h - r);
                                    pressureForce -= rijNorm * ParticleMass *
                                        (pi.Pressure + pj.Pressure) / (2 * pj.Density) * pressureGrad;

                                    // Viscosity force (Laplacian)
                                    double viscLaplacian = viscosityCoeff * (h - r);
                                    viscosityForce += (pj.Velocity - pi.Velocity) * ParticleMass *
                                        viscLaplacian / pj.Density * Viscosity;
                                }
                            }
                        }
                    }
                }

                // Gravity
                var gravityForce = Gravity * pi.Density;

                pi.Force = pressureForce + viscosityForce + gravityForce;
                _particles[i] = pi;
            }
        }

        private void ApplyExternalForces()
        {
            foreach (var force in _externalForces)
            {
                if (!force.Enabled)
                    continue;

                for (int i = 0; i < _particles.Count; i++)
                {
                    if (!_particles[i].IsActive)
                        continue;

                    var pi = _particles[i];
                    var f = force.Calculate(pi.Position, pi.Velocity, ParticleMass);
                    pi.Force += f;
                    _particles[i] = pi;
                }
            }
        }

        private void Integrate(double dt)
        {
            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive)
                    continue;

                var pi = _particles[i];

                // Acceleration
                var acceleration = pi.Force / pi.Density;

                // Semi-implicit Euler
                pi.Velocity += acceleration * dt;
                pi.Position += pi.Velocity * dt;

                _particles[i] = pi;
            }
        }

        private void HandleBoundaryCollisions()
        {
            double damping = 1.0 - BoundaryRestitution;

            for (int i = 0; i < _particles.Count; i++)
            {
                if (!_particles[i].IsActive)
                    continue;

                var pi = _particles[i];
                double radius = SmoothingRadius * 0.5;

                if (pi.Position.X - radius < Bounds.Min.X)
                {
                    pi.Position.X = Bounds.Min.X + radius;
                    pi.Velocity.X *= -BoundaryRestitution;
                }
                else if (pi.Position.X + radius > Bounds.Max.X)
                {
                    pi.Position.X = Bounds.Max.X - radius;
                    pi.Velocity.X *= -BoundaryRestitution;
                }

                if (pi.Position.Y - radius < Bounds.Min.Y)
                {
                    pi.Position.Y = Bounds.Min.Y + radius;
                    pi.Velocity.Y *= -BoundaryRestitution;
                }
                else if (pi.Position.Y + radius > Bounds.Max.Y)
                {
                    pi.Position.Y = Bounds.Max.Y - radius;
                    pi.Velocity.Y *= -BoundaryRestitution;
                }

                if (pi.Position.Z - radius < Bounds.Min.Z)
                {
                    pi.Position.Z = Bounds.Min.Z + radius;
                    pi.Velocity.Z *= -BoundaryRestitution;
                }
                else if (pi.Position.Z + radius > Bounds.Max.Z)
                {
                    pi.Position.Z = Bounds.Max.Z - radius;
                    pi.Velocity.Z *= -BoundaryRestitution;
                }

                _particles[i] = pi;
            }
        }

        private void HandleColliderCollisions()
        {
            double radius = SmoothingRadius * 0.5;

            foreach (var collider in _colliders)
            {
                if (!collider.IsActive)
                    continue;

                // Simple AABB check first
                var colliderBounds = collider.BoundingBox.Expanded(radius);

                for (int i = 0; i < _particles.Count; i++)
                {
                    if (!_particles[i].IsActive)
                        continue;

                    var pi = _particles[i];

                    if (!colliderBounds.Contains(pi.Position))
                        continue;

                    // Detailed collision check
                    if (collider is RigidBody rb)
                    {
                        HandleRigidBodyCollision(ref pi, rb, radius);
                        _particles[i] = pi;
                    }
                }
            }
        }

        private void HandleRigidBodyCollision(ref FluidParticle particle, RigidBody body, double radius)
        {
            Vector3D closestPoint;
            Vector3D normal;

            switch (body.ShapeType)
            {
                case CollisionShapeType.Sphere:
                    var toParticle = particle.Position - body.Position;
                    double dist = toParticle.Magnitude;
                    if (dist < body.Radius + radius)
                    {
                        normal = dist > PhysicsConstants.Epsilon ? toParticle / dist : Vector3D.Up;
                        closestPoint = body.Position + normal * body.Radius;

                        // Push particle out
                        particle.Position = closestPoint + normal * radius;

                        // Reflect velocity
                        double velNormal = Vector3D.Dot(particle.Velocity, normal);
                        if (velNormal < 0)
                        {
                            particle.Velocity -= normal * velNormal * (1 + CollisionRestitution);

                            // Transfer momentum to rigid body if dynamic
                            if (body.BodyType == BodyType.Dynamic)
                            {
                                body.ApplyImpulseAtPoint(
                                    normal * velNormal * ParticleMass * CollisionRestitution,
                                    closestPoint
                                );
                            }
                        }
                    }
                    break;

                case CollisionShapeType.Box:
                    // Transform to local space
                    var localPos = body.Transform.InverseTransformPoint(particle.Position);
                    var halfExtents = body.HalfExtents;

                    // Clamp to box
                    var clamped = new Vector3D(
                        Math.Clamp(localPos.X, -halfExtents.X, halfExtents.X),
                        Math.Clamp(localPos.Y, -halfExtents.Y, halfExtents.Y),
                        Math.Clamp(localPos.Z, -halfExtents.Z, halfExtents.Z)
                    );

                    var localDelta = localPos - clamped;
                    double localDist = localDelta.Magnitude;

                    if (localDist < radius || localPos == clamped) // Inside or on surface
                    {
                        if (localDist < PhysicsConstants.Epsilon)
                        {
                            // Find nearest face
                            var dists = new double[]
                            {
                                halfExtents.X - Math.Abs(localPos.X),
                                halfExtents.Y - Math.Abs(localPos.Y),
                                halfExtents.Z - Math.Abs(localPos.Z)
                            };

                            int minAxis = 0;
                            for (int a = 1; a < 3; a++)
                            {
                                if (dists[a] < dists[minAxis])
                                    minAxis = a;
                            }

                            var localNormal = Vector3D.Zero;
                            switch (minAxis)
                            {
                                case 0: localNormal = new Vector3D(Math.Sign(localPos.X), 0, 0); break;
                                case 1: localNormal = new Vector3D(0, Math.Sign(localPos.Y), 0); break;
                                case 2: localNormal = new Vector3D(0, 0, Math.Sign(localPos.Z)); break;
                            }

                            normal = body.Transform.TransformDirection(localNormal);
                            particle.Position = body.Transform.TransformPoint(clamped) + normal * radius;
                        }
                        else
                        {
                            var localNormal = localDelta / localDist;
                            normal = body.Transform.TransformDirection(localNormal);
                            closestPoint = body.Transform.TransformPoint(clamped);
                            particle.Position = closestPoint + normal * radius;
                        }

                        // Reflect velocity
                        double velNormal = Vector3D.Dot(particle.Velocity, normal);
                        if (velNormal < 0)
                        {
                            particle.Velocity -= normal * velNormal * (1 + CollisionRestitution);

                            if (body.BodyType == BodyType.Dynamic)
                            {
                                body.ApplyImpulseAtPoint(
                                    normal * velNormal * ParticleMass * CollisionRestitution,
                                    particle.Position - normal * radius
                                );
                            }
                        }
                    }
                    break;
            }
        }

        #endregion

        #region Fluent API

        /// <summary>
        /// Sets the rest density.
        /// </summary>
        public FluidSimulation WithDensity(double density)
        {
            RestDensity = density;
            return this;
        }

        /// <summary>
        /// Sets the viscosity.
        /// </summary>
        public FluidSimulation WithViscosity(double viscosity)
        {
            Viscosity = viscosity;
            return this;
        }

        /// <summary>
        /// Sets the surface tension.
        /// </summary>
        public FluidSimulation WithSurfaceTension(double tension)
        {
            SurfaceTension = tension;
            return this;
        }

        /// <summary>
        /// Sets the gas constant (stiffness).
        /// </summary>
        public FluidSimulation WithStiffness(double stiffness)
        {
            GasConstant = stiffness;
            return this;
        }

        /// <summary>
        /// Sets the gravity.
        /// </summary>
        public FluidSimulation WithGravity(Vector3D gravity)
        {
            Gravity = gravity;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a water simulation.
        /// </summary>
        public static FluidSimulation Water(AABB bounds)
            => new FluidSimulation(bounds)
                .WithDensity(1000)
                .WithViscosity(0.001)
                .WithSurfaceTension(0.0728);

        /// <summary>
        /// Creates an oil simulation.
        /// </summary>
        public static FluidSimulation Oil(AABB bounds)
            => new FluidSimulation(bounds)
                .WithDensity(900)
                .WithViscosity(0.1)
                .WithSurfaceTension(0.03);

        /// <summary>
        /// Creates a honey simulation.
        /// </summary>
        public static FluidSimulation Honey(AABB bounds)
            => new FluidSimulation(bounds)
                .WithDensity(1400)
                .WithViscosity(10.0)
                .WithSurfaceTension(0.05);

        /// <summary>
        /// Creates a lava simulation.
        /// </summary>
        public static FluidSimulation Lava(AABB bounds)
            => new FluidSimulation(bounds)
                .WithDensity(2500)
                .WithViscosity(100.0)
                .WithSurfaceTension(0.4);

        #endregion
    }

    /// <summary>
    /// Represents a fluid particle for SPH simulation.
    /// </summary>
    public struct FluidParticle
    {
        public Vector3D Position;
        public Vector3D Velocity;
        public Vector3D Force;
        public double Density;
        public double Pressure;
        public bool IsActive;
    }
}
