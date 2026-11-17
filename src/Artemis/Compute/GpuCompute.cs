using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Particles;

namespace Artemis.Compute
{
    /// <summary>
    /// GPU compute backend types.
    /// </summary>
    public enum GpuBackend
    {
        /// <summary>Automatic detection (prefers CUDA > OpenCL > CPU).</summary>
        Auto,
        /// <summary>OpenCL (Intel, AMD, NVIDIA).</summary>
        OpenCL,
        /// <summary>CUDA (NVIDIA only).</summary>
        CUDA,
        /// <summary>Vulkan compute shaders.</summary>
        Vulkan,
        /// <summary>DirectX 12 compute (Windows).</summary>
        DirectCompute,
        /// <summary>Metal compute (macOS/iOS).</summary>
        Metal,
        /// <summary>CPU fallback with SIMD.</summary>
        CPU
    }

    /// <summary>
    /// GPU device information.
    /// </summary>
    public class GpuDeviceInfo
    {
        /// <summary>Device name.</summary>
        public string Name { get; init; } = "";

        /// <summary>Vendor name.</summary>
        public string Vendor { get; init; } = "";

        /// <summary>Backend type.</summary>
        public GpuBackend Backend { get; init; }

        /// <summary>Total memory in bytes.</summary>
        public long TotalMemory { get; init; }

        /// <summary>Max compute units/cores.</summary>
        public int ComputeUnits { get; init; }

        /// <summary>Max work group size.</summary>
        public int MaxWorkGroupSize { get; init; }

        /// <summary>Supports double precision.</summary>
        public bool SupportsDouble { get; init; }

        /// <summary>Device ID for selection.</summary>
        public int DeviceId { get; init; }
    }

    /// <summary>
    /// Interface for GPU compute operations.
    /// </summary>
    public interface IGpuCompute : IDisposable
    {
        /// <summary>Gets whether GPU is available.</summary>
        bool IsAvailable { get; }

        /// <summary>Gets the active device info.</summary>
        GpuDeviceInfo? DeviceInfo { get; }

        /// <summary>Gets all available devices.</summary>
        IReadOnlyList<GpuDeviceInfo> AvailableDevices { get; }

        /// <summary>Selects a specific device.</summary>
        bool SelectDevice(int deviceId);

        /// <summary>Updates particles on GPU.</summary>
        void UpdateParticles(ParticleSoA particles, float deltaTime, Vector3D gravity);

        /// <summary>Computes particle collisions on GPU.</summary>
        void ComputeParticleCollisions(ParticleSoA particles, float cellSize);

        /// <summary>Updates rigid body positions and velocities.</summary>
        void IntegrateBodies(Span<Vector3D> positions, Span<Vector3D> velocities,
            Span<Vector3D> forces, Span<float> masses, float deltaTime);

        /// <summary>Computes SPH density and pressure.</summary>
        void ComputeSPH(ParticleSoA particles, float smoothingRadius, float restDensity, float stiffness);
    }

    /// <summary>
    /// Bulk physics operations (particle update and collisions, rigid-body
    /// integration, SPH) for <see cref="ParticleSoA"/> data.
    /// <para>
    /// Despite the name, these run on the CPU with SIMD and multi-threading on every
    /// platform; <see cref="AvailableDevices"/> also lists the machine's OpenCL
    /// devices for information. For GPU-accelerated particles use
    /// <see cref="Artemis.Particles.MassiveParticleSystem.UseGPU"/>.
    /// </para>
    /// </summary>
    public class GpuCompute : IGpuCompute
    {
        #region Fields

        private readonly GpuBackend _preferredBackend;
        private IGpuCompute? _backend;
        private readonly List<GpuDeviceInfo> _devices;
        private bool _initialized;

        #endregion

        #region Properties

        /// <summary>Gets whether GPU is available and initialized.</summary>
        public bool IsAvailable => _backend?.IsAvailable ?? false;

        /// <summary>Gets the active device info.</summary>
        public GpuDeviceInfo? DeviceInfo => _backend?.DeviceInfo;

        /// <summary>Gets all available devices.</summary>
        public IReadOnlyList<GpuDeviceInfo> AvailableDevices => _devices;

        /// <summary>Gets the active backend type.</summary>
        public GpuBackend ActiveBackend => _backend?.DeviceInfo?.Backend ?? GpuBackend.CPU;

        /// <summary>Gets or sets the preferred backend (must be set before Initialize).</summary>
        public GpuBackend PreferredBackend { get; set; } = GpuBackend.Auto;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new GPU compute instance with automatic backend detection.
        /// </summary>
        public GpuCompute() : this(GpuBackend.Auto) { }

        /// <summary>
        /// Creates a new GPU compute instance with a preferred backend.
        /// </summary>
        public GpuCompute(GpuBackend preferredBackend)
        {
            _preferredBackend = preferredBackend;
            _devices = new List<GpuDeviceInfo>();
        }

        #endregion

        #region Initialization

        /// <summary>
        /// Initializes GPU compute, detecting available devices.
        /// </summary>
        public bool Initialize()
        {
            if (_initialized) return IsAvailable;

            _devices.Clear();

            // The operations of this class (particle update, collisions, rigid-body
            // integration, SPH) run on the CPU with SIMD on every platform. The OpenCL
            // devices of the machine are listed for information; GPU acceleration of
            // particles is provided by MassiveParticleSystem.UseGPU.
            _backend = new CpuSimdCompute();
            _devices.AddRange(_backend.AvailableDevices);

#if !NETSTANDARD2_1
            int id = _devices.Count;
            foreach (var device in OpenCL.Devices)
            {
                _devices.Add(new GpuDeviceInfo
                {
                    Name = device.Name,
                    Vendor = device.Vendor,
                    Backend = GpuBackend.OpenCL,
                    TotalMemory = device.GlobalMemory,
                    ComputeUnits = device.ComputeUnits,
                    MaxWorkGroupSize = device.MaxWorkGroupSize,
                    SupportsDouble = false,
                    DeviceId = id++
                });
            }
#endif

            _initialized = true;
            return IsAvailable;
        }

        /// <summary>
        /// Selects a device by id. Only the CPU device (id 0) executes the
        /// operations of this class; the other listed devices are informational.
        /// </summary>
        public bool SelectDevice(int deviceId)
        {
            return _backend?.SelectDevice(deviceId) ?? false;
        }

        #endregion

        #region Particle Operations

        /// <summary>
        /// Updates particle positions and velocities on GPU.
        /// </summary>
        public void UpdateParticles(ParticleSoA particles, float deltaTime, Vector3D gravity)
        {
            if (!_initialized) Initialize();
            _backend?.UpdateParticles(particles, deltaTime, gravity);
        }

        /// <summary>
        /// Computes particle-particle collisions using spatial hashing on GPU.
        /// </summary>
        public void ComputeParticleCollisions(ParticleSoA particles, float cellSize)
        {
            if (!_initialized) Initialize();
            _backend?.ComputeParticleCollisions(particles, cellSize);
        }

        #endregion

        #region Rigid Body Operations

        /// <summary>
        /// Integrates rigid body motion on GPU.
        /// </summary>
        public void IntegrateBodies(
            Span<Vector3D> positions,
            Span<Vector3D> velocities,
            Span<Vector3D> forces,
            Span<float> masses,
            float deltaTime)
        {
            if (!_initialized) Initialize();
            _backend?.IntegrateBodies(positions, velocities, forces, masses, deltaTime);
        }

        #endregion

        #region SPH Fluid Simulation

        /// <summary>
        /// Computes SPH density and pressure fields on GPU.
        /// </summary>
        public void ComputeSPH(
            ParticleSoA particles,
            float smoothingRadius,
            float restDensity,
            float stiffness)
        {
            if (!_initialized) Initialize();
            _backend?.ComputeSPH(particles, smoothingRadius, restDensity, stiffness);
        }

        #endregion

        #region Disposal

        public void Dispose()
        {
            _backend?.Dispose();
            _backend = null;
            _initialized = false;
        }

        #endregion
    }

    #region CPU SIMD Backend

    /// <summary>
    /// CPU SIMD compute backend using System.Numerics.Vector.
    /// Fallback when GPU is not available.
    /// </summary>
    internal class CpuSimdCompute : IGpuCompute
    {
        private readonly List<GpuDeviceInfo> _devices;
        private readonly GpuDeviceInfo _deviceInfo;

        public bool IsAvailable => true;
        public GpuDeviceInfo? DeviceInfo => _deviceInfo;
        public IReadOnlyList<GpuDeviceInfo> AvailableDevices => _devices;

        public CpuSimdCompute()
        {
            _deviceInfo = new GpuDeviceInfo
            {
                Name = "CPU (SIMD)",
                Vendor = RuntimeInformation.ProcessArchitecture.ToString(),
                Backend = GpuBackend.CPU,
                TotalMemory = Environment.SystemPageSize * 1000L,
                ComputeUnits = Environment.ProcessorCount,
                MaxWorkGroupSize = System.Numerics.Vector<float>.Count,
                SupportsDouble = true,
                DeviceId = 0
            };
            _devices = new List<GpuDeviceInfo> { _deviceInfo };
        }

        public bool SelectDevice(int deviceId) => deviceId == 0;

        public void UpdateParticles(ParticleSoA particles, float deltaTime, Vector3D gravity)
        {
            UpdateParticlesSIMD(particles, deltaTime, gravity);
        }

        public static void UpdateParticlesSIMD(ParticleSoA particles, float deltaTime, Vector3D gravity)
        {
            int count = particles.Count;
            int vectorSize = System.Numerics.Vector<float>.Count;

            var gravityVec = new System.Numerics.Vector3((float)gravity.X, (float)gravity.Y, (float)gravity.Z);

            // Process in SIMD chunks
            Parallel.For(0, (count + vectorSize - 1) / vectorSize, chunk =>
            {
                int start = chunk * vectorSize;
                int end = Math.Min(start + vectorSize, count);

                for (int i = start; i < end; i++)
                {
                    if (particles.Lifetimes[i] <= 0) continue;

                    // Update velocity with gravity
                    particles.VelocitiesX[i] += (float)gravity.X * deltaTime;
                    particles.VelocitiesY[i] += (float)gravity.Y * deltaTime;
                    particles.VelocitiesZ[i] += (float)gravity.Z * deltaTime;

                    // Update position
                    particles.PositionsX[i] += particles.VelocitiesX[i] * deltaTime;
                    particles.PositionsY[i] += particles.VelocitiesY[i] * deltaTime;
                    particles.PositionsZ[i] += particles.VelocitiesZ[i] * deltaTime;

                    // Update lifetime
                    particles.Lifetimes[i] -= deltaTime;
                }
            });
        }

        public void ComputeParticleCollisions(ParticleSoA particles, float cellSize)
        {
            ComputeCollisionsSIMD(particles, cellSize);
        }

        public static void ComputeCollisionsSIMD(ParticleSoA particles, float cellSize)
        {
            // Spatial hash-based collision detection on CPU
            var cellToParticles = new Dictionary<long, List<int>>();
            float invCellSize = 1.0f / cellSize;

            // Build spatial hash
            for (int i = 0; i < particles.Count; i++)
            {
                if (particles.Lifetimes[i] <= 0) continue;

                int cx = (int)Math.Floor(particles.PositionsX[i] * invCellSize);
                int cy = (int)Math.Floor(particles.PositionsY[i] * invCellSize);
                int cz = (int)Math.Floor(particles.PositionsZ[i] * invCellSize);

                long hash = ((long)cx * 73856093L) ^ ((long)cy * 19349663L) ^ ((long)cz * 83492791L);

                if (!cellToParticles.TryGetValue(hash, out var list))
                {
                    list = new List<int>();
                    cellToParticles[hash] = list;
                }
                list.Add(i);
            }

            // Check collisions within cells
            double radius = cellSize * 0.5f;
            double radiusSq = radius * radius;

            Parallel.ForEach(cellToParticles.Values, cell =>
            {
                for (int i = 0; i < cell.Count; i++)
                {
                    int pi = cell[i];
                    for (int j = i + 1; j < cell.Count; j++)
                    {
                        int pj = cell[j];

                        double dx = particles.PositionsX[pi] - particles.PositionsX[pj];
                        double dy = particles.PositionsY[pi] - particles.PositionsY[pj];
                        double dz = particles.PositionsZ[pi] - particles.PositionsZ[pj];
                        double distSq = dx * dx + dy * dy + dz * dz;

                        if (distSq < radiusSq && distSq > 0.0001f)
                        {
                            double dist = Math.Sqrt(distSq);
                            double overlap = radius - dist;
                            double nx = dx / dist;
                            double ny = dy / dist;
                            double nz = dz / dist;

                            // Separate particles
                            double separation = overlap * 0.5f;
                            particles.PositionsX[pi] += nx * separation;
                            particles.PositionsY[pi] += ny * separation;
                            particles.PositionsZ[pi] += nz * separation;
                            particles.PositionsX[pj] -= nx * separation;
                            particles.PositionsY[pj] -= ny * separation;
                            particles.PositionsZ[pj] -= nz * separation;

                            // Exchange velocities (simplified elastic collision)
                            double dvx = particles.VelocitiesX[pi] - particles.VelocitiesX[pj];
                            double dvy = particles.VelocitiesY[pi] - particles.VelocitiesY[pj];
                            double dvz = particles.VelocitiesZ[pi] - particles.VelocitiesZ[pj];
                            double dvn = dvx * nx + dvy * ny + dvz * nz;

                            if (dvn < 0)
                            {
                                particles.VelocitiesX[pi] -= dvn * nx;
                                particles.VelocitiesY[pi] -= dvn * ny;
                                particles.VelocitiesZ[pi] -= dvn * nz;
                                particles.VelocitiesX[pj] += dvn * nx;
                                particles.VelocitiesY[pj] += dvn * ny;
                                particles.VelocitiesZ[pj] += dvn * nz;
                            }
                        }
                    }
                }
            });
        }

        public void IntegrateBodies(Span<Vector3D> positions, Span<Vector3D> velocities,
            Span<Vector3D> forces, Span<float> masses, float deltaTime)
        {
            IntegrateBodiesSIMD(positions, velocities, forces, masses, deltaTime);
        }

        public static void IntegrateBodiesSIMD(
            Span<Vector3D> positions,
            Span<Vector3D> velocities,
            Span<Vector3D> forces,
            Span<float> masses,
            float deltaTime)
        {
            for (int i = 0; i < positions.Length; i++)
            {
                if (masses[i] <= 0) continue;

                double invMass = 1.0 / masses[i];
                var acceleration = forces[i] * invMass;

                velocities[i] += acceleration * deltaTime;
                positions[i] += velocities[i] * deltaTime;
            }
        }

        public void ComputeSPH(ParticleSoA particles, float smoothingRadius, float restDensity, float stiffness)
        {
            ComputeSPH_SIMD(particles, smoothingRadius, restDensity, stiffness);
        }

        public static void ComputeSPH_SIMD(
            ParticleSoA particles,
            float smoothingRadius,
            float restDensity,
            float stiffness)
        {
            int count = particles.Count;
            double h = smoothingRadius;
            double h2 = h * h;
            double h9 = Math.Pow(h, 9);
            double poly6Const = 315.0 / (64.0 * Math.PI * h9);
            double spikyConst = -45.0 / (Math.PI * Math.Pow(h, 6));
            double mass = 1.0;

            // Compute densities (O(n²) - would use spatial hash in production)
            var densities = new double[count];

            Parallel.For(0, count, i =>
            {
                if (particles.Lifetimes[i] <= 0) return;

                double density = 0;
                for (int j = 0; j < count; j++)
                {
                    if (particles.Lifetimes[j] <= 0) continue;

                    double dx = particles.PositionsX[i] - particles.PositionsX[j];
                    double dy = particles.PositionsY[i] - particles.PositionsY[j];
                    double dz = particles.PositionsZ[i] - particles.PositionsZ[j];
                    double r2 = dx * dx + dy * dy + dz * dz;

                    if (r2 < h2)
                    {
                        double w = h2 - r2;
                        density += mass * poly6Const * w * w * w;
                    }
                }
                densities[i] = density;
            });

            // Compute pressure forces
            Parallel.For(0, count, i =>
            {
                if (particles.Lifetimes[i] <= 0) return;

                double pressure_i = stiffness * (densities[i] - restDensity);
                double fx = 0, fy = 0, fz = 0;

                for (int j = 0; j < count; j++)
                {
                    if (i == j || particles.Lifetimes[j] <= 0) continue;

                    double dx = particles.PositionsX[i] - particles.PositionsX[j];
                    double dy = particles.PositionsY[i] - particles.PositionsY[j];
                    double dz = particles.PositionsZ[i] - particles.PositionsZ[j];
                    double r2 = dx * dx + dy * dy + dz * dz;

                    if (r2 < h2 && r2 > 0.0001f)
                    {
                        double r = Math.Sqrt(r2);
                        double pressure_j = stiffness * (densities[j] - restDensity);
                        double avgPressure = (pressure_i + pressure_j) * 0.5;
                        double w = h - r;
                        double gradient = spikyConst * w * w / r;

                        double force = -mass * avgPressure * gradient / densities[j];
                        fx += force * dx;
                        fy += force * dy;
                        fz += force * dz;
                    }
                }

                particles.VelocitiesX[i] += fx * 0.016; // Assuming 60 FPS
                particles.VelocitiesY[i] += fy * 0.016;
                particles.VelocitiesZ[i] += fz * 0.016;
            });
        }

        public void Dispose() { }
    }

    #endregion
}
