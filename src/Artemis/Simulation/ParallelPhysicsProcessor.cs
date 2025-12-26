using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Bodies;
using Artemis.Forces;
using Artemis.Particles;

namespace Artemis.Simulation
{
    /// <summary>
    /// High-performance parallel physics processor for large-scale simulations.
    /// Uses SIMD, parallel processing, and cache-friendly data layouts.
    /// </summary>
    public class ParallelPhysicsProcessor
    {
        #region Fields

        private readonly int _threadCount;
        private readonly ParallelOptions _parallelOptions;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the number of threads used.
        /// </summary>
        public int ThreadCount => _threadCount;

        /// <summary>
        /// Gets or sets the minimum particles per batch for parallelization.
        /// </summary>
        public int MinBatchSize { get; set; } = 1000;

        /// <summary>
        /// Gets or sets whether to use SIMD operations where available.
        /// </summary>
        public bool UseSIMD { get; set; } = true;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a parallel physics processor.
        /// </summary>
        /// <param name="threadCount">Number of threads (0 = auto).</param>
        public ParallelPhysicsProcessor(int threadCount = 0)
        {
            _threadCount = threadCount > 0 ? threadCount : Environment.ProcessorCount;
            _parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = _threadCount
            };
        }

        #endregion

        #region Particle Processing

        /// <summary>
        /// Updates particles in parallel using Structure of Arrays (SoA) layout.
        /// </summary>
        public void UpdateParticlesSoA(
            ParticleSoA particles,
            double deltaTime,
            Vector3D gravity,
            IReadOnlyList<IForce>? forces = null)
        {
            int count = particles.Count;

            if (count < MinBatchSize)
            {
                // Sequential for small counts
                UpdateParticleRangeSequential(particles, 0, count, deltaTime, gravity, forces);
                return;
            }

            // Parallel processing
            int batchSize = Math.Max(MinBatchSize, count / _threadCount);

            Parallel.For(0, (count + batchSize - 1) / batchSize, _parallelOptions, batch =>
            {
                int start = batch * batchSize;
                int end = Math.Min(start + batchSize, count);
                UpdateParticleRangeSequential(particles, start, end, deltaTime, gravity, forces);
            });
        }

        private void UpdateParticleRangeSequential(
            ParticleSoA particles,
            int start,
            int end,
            double deltaTime,
            Vector3D gravity,
            IReadOnlyList<IForce>? forces)
        {
            for (int i = start; i < end; i++)
            {
                if (!particles.IsAlive[i])
                    continue;

                // Apply gravity
                particles.VelX[i] += gravity.X * deltaTime;
                particles.VelY[i] += gravity.Y * deltaTime;
                particles.VelZ[i] += gravity.Z * deltaTime;

                // Apply forces
                if (forces != null)
                {
                    var pos = new Vector3D(particles.PosX[i], particles.PosY[i], particles.PosZ[i]);
                    var vel = new Vector3D(particles.VelX[i], particles.VelY[i], particles.VelZ[i]);

                    foreach (var force in forces)
                    {
                        if (!force.Enabled) continue;
                        var f = force.Calculate(pos, vel, particles.Mass[i]);
                        double invMass = 1.0 / particles.Mass[i];
                        particles.VelX[i] += f.X * invMass * deltaTime;
                        particles.VelY[i] += f.Y * invMass * deltaTime;
                        particles.VelZ[i] += f.Z * invMass * deltaTime;
                    }
                }

                // Update position
                particles.PosX[i] += particles.VelX[i] * deltaTime;
                particles.PosY[i] += particles.VelY[i] * deltaTime;
                particles.PosZ[i] += particles.VelZ[i] * deltaTime;

                // Update lifetime
                particles.Lifetime[i] -= deltaTime;
                if (particles.Lifetime[i] <= 0)
                {
                    particles.IsAlive[i] = false;
                }
            }
        }

        /// <summary>
        /// Updates particles using SIMD (if supported).
        /// </summary>
        public void UpdateParticlesSIMD(
            ParticleSoA particles,
            double deltaTime,
            Vector3D gravity)
        {
            if (!Vector.IsHardwareAccelerated || !UseSIMD)
            {
                UpdateParticlesSoA(particles, deltaTime, gravity);
                return;
            }

            int count = particles.Count;
            int vectorSize = Vector<double>.Count;
            int vectorizedEnd = count - (count % vectorSize);

            // Vectorized gravity
            var gravX = new Vector<double>(gravity.X * deltaTime);
            var gravY = new Vector<double>(gravity.Y * deltaTime);
            var gravZ = new Vector<double>(gravity.Z * deltaTime);
            var dt = new Vector<double>(deltaTime);

            // Process in SIMD batches
            Parallel.For(0, vectorizedEnd / vectorSize, _parallelOptions, batch =>
            {
                int i = batch * vectorSize;

                // Load velocities
                var velX = new Vector<double>(particles.VelX, i);
                var velY = new Vector<double>(particles.VelY, i);
                var velZ = new Vector<double>(particles.VelZ, i);

                // Apply gravity
                velX += gravX;
                velY += gravY;
                velZ += gravZ;

                // Store velocities
                velX.CopyTo(particles.VelX, i);
                velY.CopyTo(particles.VelY, i);
                velZ.CopyTo(particles.VelZ, i);

                // Load positions
                var posX = new Vector<double>(particles.PosX, i);
                var posY = new Vector<double>(particles.PosY, i);
                var posZ = new Vector<double>(particles.PosZ, i);

                // Update positions
                posX += velX * dt;
                posY += velY * dt;
                posZ += velZ * dt;

                // Store positions
                posX.CopyTo(particles.PosX, i);
                posY.CopyTo(particles.PosY, i);
                posZ.CopyTo(particles.PosZ, i);

                // Update lifetime
                var lifetime = new Vector<double>(particles.Lifetime, i);
                lifetime -= dt;
                lifetime.CopyTo(particles.Lifetime, i);
            });

            // Process remainder
            for (int i = vectorizedEnd; i < count; i++)
            {
                if (!particles.IsAlive[i]) continue;

                particles.VelX[i] += gravity.X * deltaTime;
                particles.VelY[i] += gravity.Y * deltaTime;
                particles.VelZ[i] += gravity.Z * deltaTime;

                particles.PosX[i] += particles.VelX[i] * deltaTime;
                particles.PosY[i] += particles.VelY[i] * deltaTime;
                particles.PosZ[i] += particles.VelZ[i] * deltaTime;

                particles.Lifetime[i] -= deltaTime;
            }

            // Update alive flags
            Parallel.For(0, count, _parallelOptions, i =>
            {
                if (particles.Lifetime[i] <= 0)
                    particles.IsAlive[i] = false;
            });
        }

        #endregion

        #region Collision Detection

        /// <summary>
        /// Performs broad-phase collision detection using spatial hashing.
        /// </summary>
        public List<(int, int)> BroadPhaseParallel(
            ParticleSoA particles,
            double cellSize)
        {
            var pairs = new ConcurrentBag<(int, int)>();
            var grid = new ConcurrentDictionary<long, List<int>>();

            int count = particles.Count;

            // Hash particles to grid cells
            Parallel.For(0, count, _parallelOptions, i =>
            {
                if (!particles.IsAlive[i]) return;

                long hash = SpatialHash(
                    particles.PosX[i],
                    particles.PosY[i],
                    particles.PosZ[i],
                    cellSize
                );

                grid.AddOrUpdate(hash,
                    _ => new List<int> { i },
                    (_, list) => { lock (list) list.Add(i); return list; }
                );
            });

            // Find pairs in same cells
            Parallel.ForEach(grid.Values, _parallelOptions, cell =>
            {
                lock (cell)
                {
                    for (int i = 0; i < cell.Count; i++)
                    {
                        for (int j = i + 1; j < cell.Count; j++)
                        {
                            pairs.Add((cell[i], cell[j]));
                        }
                    }
                }
            });

            return new List<(int, int)>(pairs);
        }

        private static long SpatialHash(double x, double y, double z, double cellSize)
        {
            int ix = (int)Math.Floor(x / cellSize);
            int iy = (int)Math.Floor(y / cellSize);
            int iz = (int)Math.Floor(z / cellSize);

            // Pack into 64-bit hash
            return ((long)(ix & 0x1FFFFF) << 42) |
                   ((long)(iy & 0x1FFFFF) << 21) |
                   (long)(iz & 0x1FFFFF);
        }

        /// <summary>
        /// Resolves particle collisions in parallel.
        /// </summary>
        public void ResolveCollisionsParallel(
            ParticleSoA particles,
            List<(int, int)> pairs,
            double restitution = 0.5)
        {
            // Note: For true parallel collision resolution, we need to avoid
            // concurrent modification of the same particle. This uses atomic-like updates.

            var velocityDeltas = new ConcurrentDictionary<int, (double dx, double dy, double dz, int count)>();

            Parallel.ForEach(pairs, _parallelOptions, pair =>
            {
                int i = pair.Item1;
                int j = pair.Item2;

                double dx = particles.PosX[j] - particles.PosX[i];
                double dy = particles.PosY[j] - particles.PosY[i];
                double dz = particles.PosZ[j] - particles.PosZ[i];

                double distSq = dx * dx + dy * dy + dz * dz;
                double radiusSum = particles.Radius[i] + particles.Radius[j];

                if (distSq >= radiusSum * radiusSum || distSq < 1e-10)
                    return;

                double dist = Math.Sqrt(distSq);
                double nx = dx / dist;
                double ny = dy / dist;
                double nz = dz / dist;

                double relVelX = particles.VelX[j] - particles.VelX[i];
                double relVelY = particles.VelY[j] - particles.VelY[i];
                double relVelZ = particles.VelZ[j] - particles.VelZ[i];

                double relVelN = relVelX * nx + relVelY * ny + relVelZ * nz;

                if (relVelN > 0) return;

                double invMassI = 1.0 / particles.Mass[i];
                double invMassJ = 1.0 / particles.Mass[j];
                double totalInvMass = invMassI + invMassJ;

                double impulse = -(1 + restitution) * relVelN / totalInvMass;

                double dvxI = -impulse * invMassI * nx;
                double dvyI = -impulse * invMassI * ny;
                double dvzI = -impulse * invMassI * nz;

                double dvxJ = impulse * invMassJ * nx;
                double dvyJ = impulse * invMassJ * ny;
                double dvzJ = impulse * invMassJ * nz;

                // Accumulate velocity changes
                velocityDeltas.AddOrUpdate(i,
                    (dvxI, dvyI, dvzI, 1),
                    (_, old) => (old.dx + dvxI, old.dy + dvyI, old.dz + dvzI, old.count + 1));

                velocityDeltas.AddOrUpdate(j,
                    (dvxJ, dvyJ, dvzJ, 1),
                    (_, old) => (old.dx + dvxJ, old.dy + dvyJ, old.dz + dvzJ, old.count + 1));
            });

            // Apply accumulated velocity changes
            Parallel.ForEach(velocityDeltas, _parallelOptions, kvp =>
            {
                int i = kvp.Key;
                var delta = kvp.Value;

                // Average if multiple collisions
                double factor = 1.0 / delta.count;
                particles.VelX[i] += delta.dx * factor;
                particles.VelY[i] += delta.dy * factor;
                particles.VelZ[i] += delta.dz * factor;
            });
        }

        #endregion

        #region Rigid Body Processing

        /// <summary>
        /// Updates rigid bodies in parallel.
        /// </summary>
        public void UpdateBodiesParallel(
            IList<IPhysicsBody> bodies,
            double deltaTime,
            Vector3D gravity,
            IReadOnlyList<IForce>? forces = null)
        {
            if (bodies.Count < 100)
            {
                // Sequential for small counts
                foreach (var body in bodies)
                {
                    UpdateSingleBody(body, deltaTime, gravity, forces);
                }
                return;
            }

            Parallel.ForEach(bodies, _parallelOptions, body =>
            {
                UpdateSingleBody(body, deltaTime, gravity, forces);
            });
        }

        private void UpdateSingleBody(
            IPhysicsBody body,
            double deltaTime,
            Vector3D gravity,
            IReadOnlyList<IForce>? forces)
        {
            if (body.BodyType == BodyType.Static || !body.IsActive || body.IsSleeping)
                return;

            // Apply gravity
            body.ApplyForce(gravity * body.Mass);

            // Apply forces
            if (forces != null)
            {
                foreach (var force in forces)
                {
                    if (!force.Enabled) continue;
                    var f = force.Calculate(body.Position, body.Velocity, body.Mass);
                    body.ApplyForce(f);
                }
            }

            // Integrate
            body.Integrate(deltaTime);
        }

        #endregion

        #region SPH Fluid Processing

        /// <summary>
        /// Computes SPH density and pressure in parallel.
        /// </summary>
        public void ComputeDensityPressureParallel(
            ParticleSoA particles,
            double[] densities,
            double[] pressures,
            double smoothingRadius,
            double restDensity,
            double stiffness)
        {
            int count = particles.Count;
            double h = smoothingRadius;
            double h2 = h * h;
            double poly6Coeff = 315.0 / (64.0 * Math.PI * Math.Pow(h, 9));

            Parallel.For(0, count, _parallelOptions, i =>
            {
                if (!particles.IsAlive[i])
                {
                    densities[i] = 0;
                    pressures[i] = 0;
                    return;
                }

                double density = 0;
                double xi = particles.PosX[i];
                double yi = particles.PosY[i];
                double zi = particles.PosZ[i];

                for (int j = 0; j < count; j++)
                {
                    if (!particles.IsAlive[j]) continue;

                    double dx = xi - particles.PosX[j];
                    double dy = yi - particles.PosY[j];
                    double dz = zi - particles.PosZ[j];
                    double r2 = dx * dx + dy * dy + dz * dz;

                    if (r2 < h2)
                    {
                        double diff = h2 - r2;
                        density += particles.Mass[j] * poly6Coeff * diff * diff * diff;
                    }
                }

                densities[i] = density;
                pressures[i] = stiffness * (density - restDensity);
            });
        }

        #endregion
    }

    /// <summary>
    /// Structure of Arrays (SoA) particle data for cache-efficient processing.
    /// </summary>
    public class ParticleSoA
    {
        public int Count { get; private set; }
        public int Capacity { get; private set; }

        public double[] PosX;
        public double[] PosY;
        public double[] PosZ;
        public double[] VelX;
        public double[] VelY;
        public double[] VelZ;
        public double[] Mass;
        public double[] Radius;
        public double[] Lifetime;
        public bool[] IsAlive;
        public uint[] Color;

        /// <summary>
        /// Creates a SoA particle container.
        /// </summary>
        public ParticleSoA(int capacity)
        {
            Capacity = capacity;
            Count = 0;

            PosX = new double[capacity];
            PosY = new double[capacity];
            PosZ = new double[capacity];
            VelX = new double[capacity];
            VelY = new double[capacity];
            VelZ = new double[capacity];
            Mass = new double[capacity];
            Radius = new double[capacity];
            Lifetime = new double[capacity];
            IsAlive = new bool[capacity];
            Color = new uint[capacity];
        }

        /// <summary>
        /// Adds a particle.
        /// </summary>
        public int Add(Vector3D position, Vector3D velocity, double mass, double radius, double lifetime, uint color = 0xFFFFFFFF)
        {
            if (Count >= Capacity)
                return -1;

            int i = Count++;
            PosX[i] = position.X;
            PosY[i] = position.Y;
            PosZ[i] = position.Z;
            VelX[i] = velocity.X;
            VelY[i] = velocity.Y;
            VelZ[i] = velocity.Z;
            Mass[i] = mass;
            Radius[i] = radius;
            Lifetime[i] = lifetime;
            IsAlive[i] = true;
            Color[i] = color;

            return i;
        }

        /// <summary>
        /// Compacts the array by removing dead particles.
        /// </summary>
        public void Compact()
        {
            int writeIndex = 0;

            for (int readIndex = 0; readIndex < Count; readIndex++)
            {
                if (IsAlive[readIndex])
                {
                    if (writeIndex != readIndex)
                    {
                        PosX[writeIndex] = PosX[readIndex];
                        PosY[writeIndex] = PosY[readIndex];
                        PosZ[writeIndex] = PosZ[readIndex];
                        VelX[writeIndex] = VelX[readIndex];
                        VelY[writeIndex] = VelY[readIndex];
                        VelZ[writeIndex] = VelZ[readIndex];
                        Mass[writeIndex] = Mass[readIndex];
                        Radius[writeIndex] = Radius[readIndex];
                        Lifetime[writeIndex] = Lifetime[readIndex];
                        IsAlive[writeIndex] = true;
                        Color[writeIndex] = Color[readIndex];
                    }
                    writeIndex++;
                }
            }

            // Mark rest as dead
            for (int i = writeIndex; i < Count; i++)
            {
                IsAlive[i] = false;
            }

            Count = writeIndex;
        }

        /// <summary>
        /// Creates from a standard ParticleSystem.
        /// </summary>
        public static ParticleSoA FromParticleSystem(ParticleSystem system)
        {
            var soa = new ParticleSoA(system.MaxParticles);

            for (int i = 0; i < system.MaxParticles; i++)
            {
                var p = system.Particles[i];
                if (p.IsAlive)
                {
                    soa.Add(p.Position, p.Velocity, p.Mass, p.Radius, p.Lifetime, p.Color);
                }
            }

            return soa;
        }

        /// <summary>
        /// Gets the alive count.
        /// </summary>
        public int AliveCount
        {
            get
            {
                int count = 0;
                for (int i = 0; i < Count; i++)
                {
                    if (IsAlive[i]) count++;
                }
                return count;
            }
        }
    }
}
