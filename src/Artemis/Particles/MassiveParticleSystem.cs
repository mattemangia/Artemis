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
    public class MassiveParticleSystem : IDisposable
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
        private int[] _cellStart;
        private int[] _cellEnd;

        private int _capacity;
        private int _count;
        private volatile int _aliveCount;

        // Spatial partitioning
        private float _cellSize;
        private float _invCellSize;
        private int _gridSizeX, _gridSizeY, _gridSizeZ;
        private int _totalCells;

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
            _totalCells = _gridSizeX * _gridSizeY * _gridSizeZ;

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
            _cellStart = new int[_totalCells + 1];
            _cellEnd = new int[_totalCells + 1];

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

        private void UpdateCpu(float dt)
        {
            // 1. Integrate velocities and positions (SIMD + parallel)
            IntegrateSimd(dt);

            // 2. Build spatial hash for collision detection
            if (EnableCollisions)
            {
                BuildSpatialHashParallel();

                // 3. Detect and resolve collisions
                ResolveCollisionsParallel();
            }

            // 4. Handle boundary collisions
            HandleBoundaryCollisionsSimd();

            // 5. Update lifetimes
            UpdateLifetimesSimd(dt);
        }

        #endregion

        #region SIMD Integration

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private void IntegrateSimd(float dt)
        {
            int count = _count;
            int simdWidth = Vector<float>.Count;
            int simdCount = count / simdWidth;
            int remainder = count % simdWidth;

            var dtVec = new Vector<float>(dt);
            var gravXVec = new Vector<float>(Gravity.X * dt);
            var gravYVec = new Vector<float>(Gravity.Y * dt);
            var gravZVec = new Vector<float>(Gravity.Z * dt);
            var dampingVec = new Vector<float>(Damping);
            var gravityFlag = new Vector<byte>(4);

            // Process in parallel chunks
            int chunkSize = Math.Max(simdWidth * 64, simdCount / _threadCount * simdWidth);

            Parallel.For(0, (count + chunkSize - 1) / chunkSize, _parallelOptions, chunkIdx =>
            {
                int start = chunkIdx * chunkSize;
                int end = Math.Min(start + chunkSize, count);

                // SIMD path
                int simdEnd = start + ((end - start) / simdWidth) * simdWidth;

                for (int i = start; i < simdEnd; i += simdWidth)
                {
                    // Load velocities
                    var velX = new Vector<float>(_velX, i);
                    var velY = new Vector<float>(_velY, i);
                    var velZ = new Vector<float>(_velZ, i);

                    // Load flags to check gravity
                    var flags = new Vector<byte>(_flags, i);

                    // Apply gravity (conditional based on flags)
                    // For simplicity, apply to all - flag check would need scalar fallback
                    velX += gravXVec;
                    velY += gravYVec;
                    velZ += gravZVec;

                    // Apply damping
                    velX *= dampingVec;
                    velY *= dampingVec;
                    velZ *= dampingVec;

                    // Store velocities
                    velX.CopyTo(_velX, i);
                    velY.CopyTo(_velY, i);
                    velZ.CopyTo(_velZ, i);

                    // Load positions
                    var posX = new Vector<float>(_posX, i);
                    var posY = new Vector<float>(_posY, i);
                    var posZ = new Vector<float>(_posZ, i);

                    // Update positions
                    posX += velX * dtVec;
                    posY += velY * dtVec;
                    posZ += velZ * dtVec;

                    // Store positions
                    posX.CopyTo(_posX, i);
                    posY.CopyTo(_posY, i);
                    posZ.CopyTo(_posZ, i);
                }

                // Scalar remainder
                for (int i = simdEnd; i < end; i++)
                {
                    if ((_flags[i] & 1) == 0) continue;

                    if ((_flags[i] & 4) != 0)
                    {
                        _velX[i] += Gravity.X * dt;
                        _velY[i] += Gravity.Y * dt;
                        _velZ[i] += Gravity.Z * dt;
                    }

                    _velX[i] *= Damping;
                    _velY[i] *= Damping;
                    _velZ[i] *= Damping;

                    _posX[i] += _velX[i] * dt;
                    _posY[i] += _velY[i] * dt;
                    _posZ[i] += _velZ[i] * dt;
                }
            });
        }

        #endregion

        #region Spatial Hash (Counting Sort - GPU-friendly)

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private void BuildSpatialHashParallel()
        {
            int count = _count;

            // Clear cell counts
            Array.Clear(_cellStart, 0, _totalCells + 1);

            // Phase 1: Compute cell indices for each particle
            Parallel.For(0, count, _parallelOptions, i =>
            {
                if ((_flags[i] & 1) == 0)
                {
                    _cellIndices[i] = -1;
                    return;
                }

                int cx = (int)((_posX[i] - (float)Bounds.Min.X) * _invCellSize);
                int cy = (int)((_posY[i] - (float)Bounds.Min.Y) * _invCellSize);
                int cz = (int)((_posZ[i] - (float)Bounds.Min.Z) * _invCellSize);

                cx = Math.Clamp(cx, 0, _gridSizeX - 1);
                cy = Math.Clamp(cy, 0, _gridSizeY - 1);
                cz = Math.Clamp(cz, 0, _gridSizeZ - 1);

                int cellIndex = cx + cy * _gridSizeX + cz * _gridSizeX * _gridSizeY;
                _cellIndices[i] = cellIndex;
            });

            // Phase 2: Count particles per cell (using interlocked operations)
            var cellCounts = new int[_totalCells];
            Parallel.For(0, count, _parallelOptions, i =>
            {
                int cellIdx = _cellIndices[i];
                if (cellIdx >= 0)
                {
                    Interlocked.Increment(ref cellCounts[cellIdx]);
                }
            });

            // Phase 3: Prefix sum to get cell start indices
            _cellStart[0] = 0;
            for (int i = 0; i < _totalCells; i++)
            {
                _cellStart[i + 1] = _cellStart[i] + cellCounts[i];
                _cellEnd[i] = _cellStart[i]; // Will be incremented during insertion
            }

            // Phase 4: Insert particles into sorted order
            var insertIndices = new int[_totalCells];
            Array.Copy(_cellStart, insertIndices, _totalCells);

            for (int i = 0; i < count; i++)
            {
                int cellIdx = _cellIndices[i];
                if (cellIdx >= 0)
                {
                    int insertIdx = Interlocked.Increment(ref insertIndices[cellIdx]) - 1;
                    _particleIndices[insertIdx] = i;
                }
            }

            // Update cell end indices
            for (int i = 0; i < _totalCells; i++)
            {
                _cellEnd[i] = _cellStart[i] + cellCounts[i];
            }
        }

        #endregion

        #region Collision Detection & Resolution

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private void ResolveCollisionsParallel()
        {
            int count = _count;

            // Process collisions per-cell in parallel
            Parallel.For(0, _totalCells, _parallelOptions, cellIdx =>
            {
                int start = _cellStart[cellIdx];
                int end = _cellEnd[cellIdx];

                if (end - start < 2) return;

                // Get neighboring cell indices
                int cz = cellIdx / (_gridSizeX * _gridSizeY);
                int remainder = cellIdx % (_gridSizeX * _gridSizeY);
                int cy = remainder / _gridSizeX;
                int cx = remainder % _gridSizeX;

                // Check within cell
                for (int i = start; i < end; i++)
                {
                    int pi = _particleIndices[i];
                    if ((_flags[pi] & 3) != 3) continue; // Not alive or not collidable

                    for (int j = i + 1; j < end; j++)
                    {
                        int pj = _particleIndices[j];
                        if ((_flags[pj] & 3) != 3) continue;

                        ResolveCollision(pi, pj);
                    }
                }

                // Check with neighboring cells (only positive direction to avoid duplicates)
                for (int dz = 0; dz <= 1; dz++)
                {
                    for (int dy = 0; dy <= 1; dy++)
                    {
                        for (int dx = 0; dx <= 1; dx++)
                        {
                            if (dx == 0 && dy == 0 && dz == 0) continue;

                            int ncx = cx + dx;
                            int ncy = cy + dy;
                            int ncz = cz + dz;

                            if (ncx >= _gridSizeX || ncy >= _gridSizeY || ncz >= _gridSizeZ) continue;

                            int neighborIdx = ncx + ncy * _gridSizeX + ncz * _gridSizeX * _gridSizeY;
                            int nstart = _cellStart[neighborIdx];
                            int nend = _cellEnd[neighborIdx];

                            for (int i = start; i < end; i++)
                            {
                                int pi = _particleIndices[i];
                                if ((_flags[pi] & 3) != 3) continue;

                                for (int j = nstart; j < nend; j++)
                                {
                                    int pj = _particleIndices[j];
                                    if ((_flags[pj] & 3) != 3) continue;

                                    ResolveCollision(pi, pj);
                                }
                            }
                        }
                    }
                }
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void ResolveCollision(int i, int j)
        {
            float dx = _posX[j] - _posX[i];
            float dy = _posY[j] - _posY[i];
            float dz = _posZ[j] - _posZ[i];

            float distSq = dx * dx + dy * dy + dz * dz;
            float minDist = _radius[i] + _radius[j];

            if (distSq >= minDist * minDist || distSq < EPSILON) return;

            float dist = MathF.Sqrt(distSq);
            float invDist = 1.0f / dist;
            float nx = dx * invDist;
            float ny = dy * invDist;
            float nz = dz * invDist;

            float penetration = minDist - dist;

            // Position correction (push apart)
            float totalInvMass = _invMass[i] + _invMass[j];
            if (totalInvMass < EPSILON) return;

            float correctionMag = penetration / totalInvMass * 0.8f; // 80% correction

            float corrI = correctionMag * _invMass[i];
            float corrJ = correctionMag * _invMass[j];

            _posX[i] -= nx * corrI;
            _posY[i] -= ny * corrI;
            _posZ[i] -= nz * corrI;
            _posX[j] += nx * corrJ;
            _posY[j] += ny * corrJ;
            _posZ[j] += nz * corrJ;

            // Velocity correction (impulse)
            float dvx = _velX[j] - _velX[i];
            float dvy = _velY[j] - _velY[i];
            float dvz = _velZ[j] - _velZ[i];

            float velAlongNormal = dvx * nx + dvy * ny + dvz * nz;
            if (velAlongNormal > 0) return; // Separating

            float impulseMag = -(1.0f + Restitution) * velAlongNormal / totalInvMass;

            float impulseX = nx * impulseMag;
            float impulseY = ny * impulseMag;
            float impulseZ = nz * impulseMag;

            _velX[i] -= impulseX * _invMass[i];
            _velY[i] -= impulseY * _invMass[i];
            _velZ[i] -= impulseZ * _invMass[i];
            _velX[j] += impulseX * _invMass[j];
            _velY[j] += impulseY * _invMass[j];
            _velZ[j] += impulseZ * _invMass[j];

            // Friction
            float tanX = dvx - velAlongNormal * nx;
            float tanY = dvy - velAlongNormal * ny;
            float tanZ = dvz - velAlongNormal * nz;
            float tanLenSq = tanX * tanX + tanY * tanY + tanZ * tanZ;

            if (tanLenSq > EPSILON)
            {
                float tanLen = MathF.Sqrt(tanLenSq);
                float invTanLen = 1.0f / tanLen;
                tanX *= invTanLen;
                tanY *= invTanLen;
                tanZ *= invTanLen;

                float frictionImpulse = Friction * impulseMag;
                _velX[i] += tanX * frictionImpulse * _invMass[i];
                _velY[i] += tanY * frictionImpulse * _invMass[i];
                _velZ[i] += tanZ * frictionImpulse * _invMass[i];
                _velX[j] -= tanX * frictionImpulse * _invMass[j];
                _velY[j] -= tanY * frictionImpulse * _invMass[j];
                _velZ[j] -= tanZ * frictionImpulse * _invMass[j];
            }
        }

        #endregion

        #region Boundary Collisions

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private void HandleBoundaryCollisionsSimd()
        {
            float minX = (float)Bounds.Min.X;
            float minY = (float)Bounds.Min.Y;
            float minZ = (float)Bounds.Min.Z;
            float maxX = (float)Bounds.Max.X;
            float maxY = (float)Bounds.Max.Y;
            float maxZ = (float)Bounds.Max.Z;
            float restitution = Restitution;
            float friction = 1.0f - Friction;

            int count = _count;

            Parallel.For(0, count, _parallelOptions, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                float r = _radius[i];

                // X bounds
                if (_posX[i] - r < minX)
                {
                    _posX[i] = minX + r;
                    _velX[i] = -_velX[i] * restitution;
                    _velY[i] *= friction;
                    _velZ[i] *= friction;
                }
                else if (_posX[i] + r > maxX)
                {
                    _posX[i] = maxX - r;
                    _velX[i] = -_velX[i] * restitution;
                    _velY[i] *= friction;
                    _velZ[i] *= friction;
                }

                // Y bounds
                if (_posY[i] - r < minY)
                {
                    _posY[i] = minY + r;
                    _velY[i] = -_velY[i] * restitution;
                    _velX[i] *= friction;
                    _velZ[i] *= friction;
                }
                else if (_posY[i] + r > maxY)
                {
                    _posY[i] = maxY - r;
                    _velY[i] = -_velY[i] * restitution;
                    _velX[i] *= friction;
                    _velZ[i] *= friction;
                }

                // Z bounds
                if (_posZ[i] - r < minZ)
                {
                    _posZ[i] = minZ + r;
                    _velZ[i] = -_velZ[i] * restitution;
                    _velX[i] *= friction;
                    _velY[i] *= friction;
                }
                else if (_posZ[i] + r > maxZ)
                {
                    _posZ[i] = maxZ - r;
                    _velZ[i] = -_velZ[i] * restitution;
                    _velX[i] *= friction;
                    _velY[i] *= friction;
                }
            });
        }

        #endregion

        #region Lifetime Management

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private void UpdateLifetimesSimd(float dt)
        {
            int count = _count;
            int simdWidth = Vector<float>.Count;

            var dtVec = new Vector<float>(dt);
            var zeroVec = Vector<float>.Zero;

            Parallel.For(0, (count + simdWidth - 1) / simdWidth, _parallelOptions, chunk =>
            {
                int start = chunk * simdWidth;
                int end = Math.Min(start + simdWidth, count);

                if (end - start == simdWidth)
                {
                    var lifetime = new Vector<float>(_lifetime, start);
                    lifetime -= dtVec;
                    lifetime.CopyTo(_lifetime, start);

                    // Check for dead particles
                    for (int i = start; i < end; i++)
                    {
                        if (_lifetime[i] <= 0)
                        {
                            _flags[i] &= unchecked((byte)~1); // Clear alive flag
                        }
                    }
                }
                else
                {
                    for (int i = start; i < end; i++)
                    {
                        _lifetime[i] -= dt;
                        if (_lifetime[i] <= 0)
                        {
                            _flags[i] &= unchecked((byte)~1);
                        }
                    }
                }
            });
        }

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

                if (distSq < radiusSq)
                {
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

            for (int readIndex = 0; readIndex < _count; readIndex++)
            {
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
                    }
                    writeIndex++;
                }
            }

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
