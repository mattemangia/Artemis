using System;
using System.Buffers;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Compute;

namespace Artemis.Particles
{
    /// <summary>
    /// Ultra high-performance particle system designed for 1M+ particles.
    /// Uses SIMD, multi-threading, and GPU acceleration (OpenCL/Metal/CUDA).
    /// Cross-platform: Windows, Linux, macOS (including Apple Silicon M1/M2/M3).
    /// </summary>
    public partial class MassiveParticleSystem : IDisposable
    {
        #region Constants

        private const int SIMD_ALIGNMENT = 64; // Cache line alignment
        private const int DEFAULT_CELL_CAPACITY = 64;
        private const float EPSILON = 1e-6f;

        #endregion

        #region Fields - Structure of Arrays (SoA) for SIMD efficiency

        // Position arrays (aligned for SIMD)
        private float[] _posX;
        private float[] _posY;
        private float[] _posZ;

        // Velocity arrays
        private float[] _velX;
        private float[] _velY;
        private float[] _velZ;

        // Previous position (for Verlet integration)
        private float[] _prevX;
        private float[] _prevY;
        private float[] _prevZ;

        // Particle properties
        private float[] _mass;
        private float[] _invMass;
        private float[] _radius;
        private float[] _lifetime;
        private uint[] _color;
        private byte[] _flags; // bit 0: alive, bit 1: collides, bit 2: affected by gravity

        // Spatial hash data
        private int[] _cellIndices;
        private int[] _particleIndices;

        private int _capacity;
        private int _count;
        private volatile int _aliveCount;

        // Spatial partitioning
        private float _cellSize;
        private float _invCellSize;
        private int _gridSizeX, _gridSizeY, _gridSizeZ;
        private long _totalCells;

        // GPU compute
        private MassiveGpuCompute? _gpuCompute;
        private bool _useGpu;
        private bool _gpuInitialized;

        // Thread management
        private readonly int _threadCount;
        private readonly ParallelOptions _parallelOptions;

        // Object pool for temporary arrays
        private readonly ArrayPool<int> _intPool;
        private readonly ArrayPool<float> _floatPool;

        #endregion

        #region Properties

        /// <summary>
        /// Maximum number of particles this system can hold.
        /// </summary>
        public int Capacity => _capacity;

        /// <summary>
        /// Current number of particles (including dead ones).
        /// </summary>
        public int Count => _count;

        /// <summary>
        /// Number of alive particles.
        /// </summary>
        public int AliveCount => _aliveCount;

        /// <summary>
        /// Gravity vector applied to particles.
        /// </summary>
        public Vector3 Gravity { get; set; } = new Vector3(0, -9.81f, 0);

        /// <summary>
        /// Simulation bounds.
        /// </summary>
        public AABB Bounds { get; set; }

        /// <summary>
        /// Collision restitution (bounciness).
        /// </summary>
        public float Restitution { get; set; } = 0.3f;

        /// <summary>
        /// Friction coefficient.
        /// </summary>
        public float Friction { get; set; } = 0.1f;

        /// <summary>
        /// Damping factor (0-1, applied each frame).
        /// </summary>
        public float Damping { get; set; } = 0.999f;

        /// <summary>
        /// Whether to use GPU acceleration.
        /// </summary>
        public bool UseGPU
        {
            get => _useGpu;
            set
            {
                if (value && !_gpuInitialized)
                {
                    InitializeGpu();
                }
                _useGpu = value && _gpuInitialized;
            }
        }

        /// <summary>
        /// Whether to enable particle-particle collisions.
        /// </summary>
        public bool EnableCollisions { get; set; } = true;

        /// <summary>
        /// Whether to use Verlet integration (more stable for constraints).
        /// </summary>
        public bool UseVerletIntegration { get; set; } = false;

        /// <summary>
        /// Gets the GPU compute backend info.
        /// </summary>
        public string GpuBackendInfo => _gpuCompute?.BackendInfo ?? "CPU Only";

        /// <summary>
        /// Direct access to position X array for rendering.
        /// </summary>
        public ReadOnlySpan<float> PositionsX => new ReadOnlySpan<float>(_posX, 0, _count);

        /// <summary>
        /// Direct access to position Y array for rendering.
        /// </summary>
        public ReadOnlySpan<float> PositionsY => new ReadOnlySpan<float>(_posY, 0, _count);

        /// <summary>
        /// Direct access to position Z array for rendering.
        /// </summary>
        public ReadOnlySpan<float> PositionsZ => new ReadOnlySpan<float>(_posZ, 0, _count);

        /// <summary>
        /// Direct access to color array for rendering.
        /// </summary>
        public ReadOnlySpan<uint> Colors => new ReadOnlySpan<uint>(_color, 0, _count);

        /// <summary>
        /// Direct access to radius array for rendering.
        /// </summary>
        public ReadOnlySpan<float> Radii => new ReadOnlySpan<float>(_radius, 0, _count);

        /// <summary>
        /// Direct access to alive flags for rendering.
        /// </summary>
        public ReadOnlySpan<byte> Flags => new ReadOnlySpan<byte>(_flags, 0, _count);

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a massive particle system optimized for large-scale simulations.
        /// </summary>
        /// <param name="capacity">Maximum number of particles (default: 1 million).</param>
        /// <param name="bounds">Simulation bounds.</param>
        /// <param name="cellSize">Spatial hash cell size (should be >= max particle diameter).</param>
        public MassiveParticleSystem(int capacity = 1_000_000, AABB? bounds = null, float cellSize = 1.0f)
        {
            _capacity = capacity;
            _count = 0;
            _aliveCount = 0;

            // Set bounds
            Bounds = bounds ?? new AABB(
                new Vector3D(-500, -500, -500),
                new Vector3D(500, 500, 500)
            );

            // Initialize spatial hash parameters
            _cellSize = cellSize;
            _invCellSize = 1.0f / cellSize;

            var boundsSize = Bounds.Max - Bounds.Min;
            _gridSizeX = Math.Max(1, (int)Math.Ceiling(boundsSize.X / cellSize));
            _gridSizeY = Math.Max(1, (int)Math.Ceiling(boundsSize.Y / cellSize));
            _gridSizeZ = Math.Max(1, (int)Math.Ceiling(boundsSize.Z / cellSize));
            _totalCells = (long)_gridSizeX * _gridSizeY * _gridSizeZ;

            // Thread configuration
            _threadCount = Environment.ProcessorCount;
            _parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = _threadCount
            };

            // Object pools
            _intPool = ArrayPool<int>.Shared;
            _floatPool = ArrayPool<float>.Shared;

            // Allocate SoA arrays
            AllocateArrays(capacity);
        }

        private void AllocateArrays(int capacity)
        {
            // Ensure capacity is aligned to SIMD width
            int simdWidth = Vector<float>.Count;
            int alignedCapacity = ((capacity + simdWidth - 1) / simdWidth) * simdWidth;

            _posX = new float[alignedCapacity];
            _posY = new float[alignedCapacity];
            _posZ = new float[alignedCapacity];
            _velX = new float[alignedCapacity];
            _velY = new float[alignedCapacity];
            _velZ = new float[alignedCapacity];
            _prevX = new float[alignedCapacity];
            _prevY = new float[alignedCapacity];
            _prevZ = new float[alignedCapacity];
            _mass = new float[alignedCapacity];
            _invMass = new float[alignedCapacity];
            _radius = new float[alignedCapacity];
            _lifetime = new float[alignedCapacity];
            _color = new uint[alignedCapacity];
            _flags = new byte[alignedCapacity];

            // Spatial hash arrays
            _cellIndices = new int[alignedCapacity];
            _particleIndices = new int[alignedCapacity];
            AllocateSolverArrays(alignedCapacity);

            _capacity = alignedCapacity;
        }

        #endregion

        #region GPU Initialization

        private void InitializeGpu()
        {
            try
            {
                _gpuCompute = new MassiveGpuCompute();
                _gpuInitialized = _gpuCompute.Initialize();

                if (_gpuInitialized)
                {
                    _gpuCompute.AllocateBuffers(_capacity);
                }
            }
            catch (Exception)
            {
                _gpuCompute = null;
                _gpuInitialized = false;
            }
        }

        #endregion

        #region Particle Management

        /// <summary>
        /// Adds a single particle.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int AddParticle(
            Vector3 position,
            Vector3 velocity,
            float mass = 1.0f,
            float radius = 0.1f,
            float lifetime = float.MaxValue,
            uint color = 0xFFFFFFFF,
            bool affectedByGravity = true,
            bool collides = true)
        {
            int index = Interlocked.Increment(ref _count) - 1;
            if (index >= _capacity)
            {
                Interlocked.Decrement(ref _count);
                return -1;
            }

            _posX[index] = position.X;
            _posY[index] = position.Y;
            _posZ[index] = position.Z;
            _velX[index] = velocity.X;
            _velY[index] = velocity.Y;
            _velZ[index] = velocity.Z;
            _prevX[index] = position.X - velocity.X * 0.016f;
            _prevY[index] = position.Y - velocity.Y * 0.016f;
            _prevZ[index] = position.Z - velocity.Z * 0.016f;
            _mass[index] = mass;
            _invMass[index] = mass > 0 ? 1.0f / mass : 0;
            _radius[index] = radius;
            _lifetime[index] = lifetime;
            _color[index] = color;

            byte flags = 1; // alive
            if (collides) flags |= 2;
            if (affectedByGravity) flags |= 4;
            _flags[index] = flags;
            _rest[index] = 0;
            _wake[index] = 0;
            if (radius > _maxRadius) _maxRadius = radius;

            Interlocked.Increment(ref _aliveCount);
            return index;
        }

        /// <summary>
        /// Spawns particles in a sphere shape. Highly optimized for bulk spawning.
        /// </summary>
        public int SpawnSphere(
            Vector3 center,
            float sphereRadius,
            float particleSpacing,
            float particleMass = 1.0f,
            float particleRadius = 0.1f,
            float lifetime = float.MaxValue,
            uint color = 0xFFFFFFFF)
        {
            // Pre-calculate positions
            var positions = new List<Vector3>();
            float spacing = particleSpacing;

            for (float x = -sphereRadius; x <= sphereRadius; x += spacing)
            {
                for (float y = -sphereRadius; y <= sphereRadius; y += spacing)
                {
                    for (float z = -sphereRadius; z <= sphereRadius; z += spacing)
                    {
                        if (x * x + y * y + z * z <= sphereRadius * sphereRadius)
                        {
                            positions.Add(new Vector3(center.X + x, center.Y + y, center.Z + z));
                        }
                    }
                }
            }

            return SpawnBatch(positions, particleMass, particleRadius, lifetime, color);
        }

        /// <summary>
        /// Spawns particles in a box shape.
        /// </summary>
        public int SpawnBox(
            Vector3 min,
            Vector3 max,
            float particleSpacing,
            float particleMass = 1.0f,
            float particleRadius = 0.1f,
            float lifetime = float.MaxValue,
            uint color = 0xFFFFFFFF)
        {
            var positions = new List<Vector3>();

            for (float x = min.X; x <= max.X; x += particleSpacing)
            {
                for (float y = min.Y; y <= max.Y; y += particleSpacing)
                {
                    for (float z = min.Z; z <= max.Z; z += particleSpacing)
                    {
                        positions.Add(new Vector3(x, y, z));
                    }
                }
            }

            return SpawnBatch(positions, particleMass, particleRadius, lifetime, color);
        }

        /// <summary>
        /// Batch spawn particles (most efficient).
        /// </summary>
        public int SpawnBatch(
            IList<Vector3> positions,
            float mass = 1.0f,
            float radius = 0.1f,
            float lifetime = float.MaxValue,
            uint color = 0xFFFFFFFF)
        {
            int count = positions.Count;
            int startIndex = Interlocked.Add(ref _count, count) - count;

            if (startIndex + count > _capacity)
            {
                Interlocked.Add(ref _count, -count);
                return 0;
            }

            float invMass = mass > 0 ? 1.0f / mass : 0;
            byte flags = (byte)(1 | 2 | 4); // alive, collides, affected by gravity
            if (radius > _maxRadius) _maxRadius = radius;

            // Parallel batch initialization
            Parallel.For(0, count, _parallelOptions, i =>
            {
                int idx = startIndex + i;
                var pos = positions[i];

                _posX[idx] = pos.X;
                _posY[idx] = pos.Y;
                _posZ[idx] = pos.Z;
                _velX[idx] = 0;
                _velY[idx] = 0;
                _velZ[idx] = 0;
                _prevX[idx] = pos.X;
                _prevY[idx] = pos.Y;
                _prevZ[idx] = pos.Z;
                _mass[idx] = mass;
                _invMass[idx] = invMass;
                _radius[idx] = radius;
                _lifetime[idx] = lifetime;
                _color[idx] = color;
                _flags[idx] = flags;
                _rest[idx] = 0;
                _wake[idx] = 0;
            });

            Interlocked.Add(ref _aliveCount, count);
            return count;
        }

        /// <summary>
        /// Clears all particles.
        /// </summary>
        public void Clear()
        {
            Array.Clear(_flags, 0, _count);
            _count = 0;
            _aliveCount = 0;
            _maxRadius = 0;
            ClearGroups();
        }

        #endregion

        #region Simulation Update

        /// <summary>
        /// Updates the particle simulation.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        /// <param name="subSteps">Number of sub-steps for stability.</param>
        public void Update(float deltaTime, int subSteps = 4)
        {
            if (_count == 0) return;

            float dt = deltaTime / subSteps;

            for (int step = 0; step < subSteps; step++)
            {
                if (_useGpu && _gpuCompute != null && _gpuCompute.IsAvailable)
                {
                    UpdateGpu(dt);
                }
                else
                {
                    UpdateCpu(dt);
                }
            }

            // Update alive count
            UpdateAliveCount();
        }

        private void UpdateGpu(float dt)
        {
            // Upload data to GPU
            _gpuCompute!.UploadParticleData(
                _posX, _posY, _posZ,
                _velX, _velY, _velZ,
                _mass, _radius, _flags,
                _count
            );

            // GPU kernels
            _gpuCompute.IntegrateParticles(dt, Gravity, Damping);

            if (EnableCollisions)
            {
                _gpuCompute.BuildSpatialHash(_cellSize, _gridSizeX, _gridSizeY, _gridSizeZ);
                _gpuCompute.DetectAndResolveCollisions(Restitution, Friction);
            }

            _gpuCompute.HandleBoundaryCollisions(
                (float)Bounds.Min.X, (float)Bounds.Min.Y, (float)Bounds.Min.Z,
                (float)Bounds.Max.X, (float)Bounds.Max.Y, (float)Bounds.Max.Z,
                Restitution
            );

            // Download results
            _gpuCompute.DownloadParticleData(
                _posX, _posY, _posZ,
                _velX, _velY, _velZ,
                _flags,
                _count
            );
        }

        #endregion

        #region Lifetime Management

        private void UpdateAliveCount()
        {
            int alive = 0;
            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & 1) != 0) alive++;
            }
            _aliveCount = alive;
        }

        #endregion

        #region Utility Methods

        /// <summary>
        /// Applies an explosion force at the given position.
        /// </summary>
        public void ApplyExplosion(Vector3 center, float force, float radius)
            => ApplyExplosion(center, force, radius, releaseFrozen: true);

        /// <summary>
        /// Applies an explosion impulse. Sleeping particles in range are woken and,
        /// if <paramref name="releaseFrozen"/> is set, frozen ones are released.
        /// </summary>
        public void ApplyExplosion(Vector3 center, float force, float radius, bool releaseFrozen)
        {
            float radiusSq = radius * radius;

            Parallel.For(0, _count, _parallelOptions, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                float dx = _posX[i] - center.X;
                float dy = _posY[i] - center.Y;
                float dz = _posZ[i] - center.Z;
                float distSq = dx * dx + dy * dy + dz * dz;

                if (distSq < radiusSq && distSq > EPSILON)
                {
                    if ((_flags[i] & FlagFrozen) != 0)
                    {
                        if (!releaseFrozen) return;
                        _flags[i] &= unchecked((byte)~FlagFrozen);
                    }
                    _flags[i] &= unchecked((byte)~FlagSleeping);
                    _rest[i] = 0;

                    float dist = MathF.Sqrt(distSq);
                    float falloff = 1.0f - dist / radius;
                    float impulse = force * falloff * _invMass[i];

                    _velX[i] += (dx / dist) * impulse;
                    _velY[i] += (dy / dist) * impulse;
                    _velZ[i] += (dz / dist) * impulse;
                }
            });
        }

        /// <summary>
        /// Applies a force to all particles within a radius.
        /// </summary>
        public void ApplyForceInRadius(Vector3 center, float radius, Vector3 force)
        {
            float radiusSq = radius * radius;

            Parallel.For(0, _count, _parallelOptions, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                float dx = _posX[i] - center.X;
                float dy = _posY[i] - center.Y;
                float dz = _posZ[i] - center.Z;
                float distSq = dx * dx + dy * dy + dz * dz;

                if (distSq < radiusSq && (_flags[i] & FlagFrozen) == 0)
                {
                    _flags[i] &= unchecked((byte)~FlagSleeping);
                    _rest[i] = 0;
                    _velX[i] += force.X * _invMass[i];
                    _velY[i] += force.Y * _invMass[i];
                    _velZ[i] += force.Z * _invMass[i];
                }
            });
        }

        /// <summary>
        /// Gets particle data for a specific index.
        /// </summary>
        public (Vector3 Position, Vector3 Velocity, float Mass, float Radius, uint Color, bool IsAlive) GetParticle(int index)
        {
            if (index < 0 || index >= _count)
                throw new ArgumentOutOfRangeException(nameof(index));

            return (
                new Vector3(_posX[index], _posY[index], _posZ[index]),
                new Vector3(_velX[index], _velY[index], _velZ[index]),
                _mass[index],
                _radius[index],
                _color[index],
                (_flags[index] & 1) != 0
            );
        }

        /// <summary>
        /// Compacts the particle array by removing dead particles.
        /// Call periodically to maintain performance.
        /// </summary>
        public void Compact()
        {
            int writeIndex = 0;
            var newIndex = new int[_count];

            for (int readIndex = 0; readIndex < _count; readIndex++)
            {
                newIndex[readIndex] = (_flags[readIndex] & 1) != 0 ? writeIndex : -1;

                if ((_flags[readIndex] & 1) != 0)
                {
                    if (writeIndex != readIndex)
                    {
                        _posX[writeIndex] = _posX[readIndex];
                        _posY[writeIndex] = _posY[readIndex];
                        _posZ[writeIndex] = _posZ[readIndex];
                        _velX[writeIndex] = _velX[readIndex];
                        _velY[writeIndex] = _velY[readIndex];
                        _velZ[writeIndex] = _velZ[readIndex];
                        _prevX[writeIndex] = _prevX[readIndex];
                        _prevY[writeIndex] = _prevY[readIndex];
                        _prevZ[writeIndex] = _prevZ[readIndex];
                        _mass[writeIndex] = _mass[readIndex];
                        _invMass[writeIndex] = _invMass[readIndex];
                        _radius[writeIndex] = _radius[readIndex];
                        _lifetime[writeIndex] = _lifetime[readIndex];
                        _color[writeIndex] = _color[readIndex];
                        _flags[writeIndex] = _flags[readIndex];
                        _rest[writeIndex] = _rest[readIndex];
                        _wake[writeIndex] = _wake[readIndex];
                        _groupOf[writeIndex] = _groupOf[readIndex];
                        _restX[writeIndex] = _restX[readIndex];
                        _restY[writeIndex] = _restY[readIndex];
                        _restZ[writeIndex] = _restZ[readIndex];
                    }
                    writeIndex++;
                }
            }

            for (int i = writeIndex; i < _count; i++)
                _groupOf[i] = -1;

            RemapGroups(newIndex);

            _count = writeIndex;
            _aliveCount = writeIndex;
        }

        #endregion

        #region Disposal

        public void Dispose()
        {
            _gpuCompute?.Dispose();
            _gpuCompute = null;
        }

        #endregion
    }
}
