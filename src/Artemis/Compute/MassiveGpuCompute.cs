using System;
using System.Collections.Generic;

namespace Artemis.Compute
{
    /// <summary>
    /// GPU backend of <see cref="Artemis.Particles.MassiveParticleSystem"/>.
    /// <para>
    /// It runs the particle solver as OpenCL kernels, which covers NVIDIA, AMD and
    /// Intel GPUs on Windows and Linux and (deprecated, but still shipped) macOS.
    /// The kernels are a port of the CPU solver, so a scene behaves the same on
    /// either path. Without an OpenCL driver <see cref="IsAvailable"/> is false and
    /// the particle system stays on the CPU.
    /// </para>
    /// </summary>
    public sealed class MassiveGpuCompute : IDisposable
    {
        private bool _initialized;
        private bool _available;
        private string _info = "CPU (no OpenCL device)";

        /// <summary>
        /// Whether a usable OpenCL device was found and the kernels compiled.
        /// </summary>
        public bool IsAvailable => _available;

        /// <summary>
        /// The backend in use: <see cref="GpuBackendType.OpenCL"/>, or
        /// <see cref="GpuBackendType.CPU_SIMD"/> when no device is available.
        /// </summary>
        public GpuBackendType ActiveBackend => _available ? GpuBackendType.OpenCL : GpuBackendType.CPU_SIMD;

        /// <summary>
        /// Human-readable description of the device.
        /// </summary>
        public string BackendInfo => _info;

        /// <summary>
        /// Why initialisation failed, if it did (no driver, no device, compile error).
        /// </summary>
        public string? InitializationError { get; private set; }

        /// <summary>
        /// Looks for an OpenCL device and compiles the particle kernels on it.
        /// </summary>
        /// <param name="allowCpuDevice">
        /// Also accept OpenCL CPU devices (POCL, the Intel CPU runtime). Off by default:
        /// on a CPU the built-in multi-threaded solver is the better choice.
        /// </param>
        /// <returns>Whether the GPU path can be used.</returns>
        public bool Initialize(bool allowCpuDevice = false)
        {
            if (_initialized)
                return _available;
            _initialized = true;

#if NETSTANDARD2_1
            InitializationError = "OpenCL is not supported on .NET Standard 2.1.";
            return false;
#else
            try
            {
                var device = OpenCL.SelectDevice(allowCpuDevice);
                if (device == null)
                {
                    InitializationError = OpenCL.IsAvailable
                        ? "Only OpenCL CPU devices were found (pass allowCpuDevice: true to use them)."
                        : "No OpenCL driver or device was found.";
                    return false;
                }

                Runtime = new OpenCLRuntime(device, ParticleKernels.Source, "-cl-std=CL1.2");
                Device = device;
                _available = true;
                _info = $"OpenCL: {device}";
            }
            catch (Exception exception) when (exception is OpenCLException || exception is DllNotFoundException || exception is EntryPointNotFoundException)
            {
                InitializationError = exception.Message;
                Runtime?.Dispose();
                Runtime = null;
            }

            return _available;
#endif
        }

        /// <summary>
        /// Kept for compatibility; buffers are now sized on demand by the particle system.
        /// </summary>
        public void AllocateBuffers(int capacity)
        {
        }

#if !NETSTANDARD2_1
        /// <summary>
        /// The OpenCL device in use, or <c>null</c>.
        /// </summary>
        public OpenCLDeviceInfo? Device { get; private set; }

        internal OpenCLRuntime? Runtime { get; private set; }

        // Device buffers, allocated by EnsureCapacity.
        internal int Capacity;
        internal IntPtr PX, PY, PZ, QX, QY, QZ, VX, VY, VZ;
        internal IntPtr InvMass, Radius, Lifetime, Flags, Rest, Wake, Body, Cell, SortedIndex;
        internal IntPtr DpX, DpY, DpZ, StX, StY, StZ;
        internal IntPtr SX, SY, SZ, SR, SW, SMX, SMY, SMZ, SKey, SBody;
        internal IntPtr Counters;

        internal int TableCapacity;
        internal IntPtr BucketCount, BucketStart, Cursor;
        internal readonly List<IntPtr> ScanSums = new();
        internal readonly List<int> ScanSumLengths = new();

        // Soft-body cluster data.
        internal int ClusterCapacity, SlotCapacity, SoftBodyCapacity;
        internal IntPtr ClusterStart, Members, RestX, RestY, RestZ, Rotation, GoalX, GoalY, GoalZ;
        internal IntPtr SlotStart, Slots, Alpha, Beta;
        internal int ChunkCapacity, BodyMemberCapacity, RigidCapacity;
        internal IntPtr BodyStart, BodyMembers, BodyChunkStart, ChunkBody, ChunkFrom, ChunkTo, ChunkPartial, BodyRigid;

        internal void EnsureCapacity(int capacity)
        {
            var rt = Runtime!;
            if (capacity <= Capacity)
                return;

            ReleaseParticleBuffers();
            Capacity = capacity;
            long f = capacity * 4L, b = capacity;

            PX = ParticleBuffer(f); PY = ParticleBuffer(f); PZ = ParticleBuffer(f);
            QX = ParticleBuffer(f); QY = ParticleBuffer(f); QZ = ParticleBuffer(f);
            VX = ParticleBuffer(f); VY = ParticleBuffer(f); VZ = ParticleBuffer(f);
            InvMass = ParticleBuffer(f); Radius = ParticleBuffer(f); Lifetime = ParticleBuffer(f);
            Flags = ParticleBuffer(b); Rest = ParticleBuffer(b);
            Wake = ParticleBuffer(f); Body = ParticleBuffer(f); Cell = ParticleBuffer(f); SortedIndex = ParticleBuffer(f);
            DpX = ParticleBuffer(f); DpY = ParticleBuffer(f); DpZ = ParticleBuffer(f);
            StX = ParticleBuffer(f); StY = ParticleBuffer(f); StZ = ParticleBuffer(f);
            SX = ParticleBuffer(f); SY = ParticleBuffer(f); SZ = ParticleBuffer(f);
            SR = ParticleBuffer(f); SW = ParticleBuffer(f);
            SMX = ParticleBuffer(f); SMY = ParticleBuffer(f); SMZ = ParticleBuffer(f);
            SKey = ParticleBuffer(f); SBody = ParticleBuffer(f);
            SlotStart = ParticleBuffer((capacity + 1) * 4L);
            Counters = ParticleBuffer(16);
        }

        internal void EnsureTable(int tableSize)
        {
            var rt = Runtime!;
            if (tableSize <= TableCapacity)
                return;

            rt.ReleaseBuffer(ref BucketCount);
            rt.ReleaseBuffer(ref BucketStart);
            rt.ReleaseBuffer(ref Cursor);
            for (int i = 0; i < ScanSums.Count; i++)
            {
                var buffer = ScanSums[i];
                rt.ReleaseBuffer(ref buffer);
            }
            ScanSums.Clear();
            ScanSumLengths.Clear();

            TableCapacity = tableSize;
            BucketCount = rt.CreateBuffer(tableSize * 4L);
            BucketStart = rt.CreateBuffer((tableSize + 4) * 4L);
            Cursor = rt.CreateBuffer(tableSize * 4L);

            // One block-sum buffer per level of the recursive prefix sum.
            int length = tableSize;
            while (length > 1)
            {
                length = (length + ParticleKernels.ScanBlock - 1) / ParticleKernels.ScanBlock;
                ScanSums.Add(rt.CreateBuffer(Math.Max(length, 1) * 4L));
                ScanSumLengths.Add(length);
            }
        }

        internal void EnsureSoftCapacity(int clusters, int slots, int bodies)
        {
            var rt = Runtime!;
            if (clusters > ClusterCapacity)
            {
                rt.ReleaseBuffer(ref ClusterStart);
                rt.ReleaseBuffer(ref Rotation);
                ClusterCapacity = Math.Max(clusters, ClusterCapacity * 2);
                ClusterStart = rt.CreateBuffer((ClusterCapacity + 1) * 4L);
                Rotation = rt.CreateBuffer(ClusterCapacity * 16L);
            }

            if (slots > SlotCapacity)
            {
                rt.ReleaseBuffer(ref Members); rt.ReleaseBuffer(ref Slots);
                rt.ReleaseBuffer(ref RestX); rt.ReleaseBuffer(ref RestY); rt.ReleaseBuffer(ref RestZ);
                rt.ReleaseBuffer(ref GoalX); rt.ReleaseBuffer(ref GoalY); rt.ReleaseBuffer(ref GoalZ);
                SlotCapacity = Math.Max(slots, SlotCapacity * 2);
                long s = SlotCapacity * 4L;
                Members = rt.CreateBuffer(s); Slots = rt.CreateBuffer(s);
                RestX = rt.CreateBuffer(s); RestY = rt.CreateBuffer(s); RestZ = rt.CreateBuffer(s);
                GoalX = rt.CreateBuffer(s); GoalY = rt.CreateBuffer(s); GoalZ = rt.CreateBuffer(s);
            }

            if (bodies > SoftBodyCapacity)
            {
                rt.ReleaseBuffer(ref Alpha);
                rt.ReleaseBuffer(ref Beta);
                SoftBodyCapacity = Math.Max(bodies, SoftBodyCapacity * 2);
                Alpha = rt.CreateBuffer(SoftBodyCapacity * 4L);
                Beta = rt.CreateBuffer(SoftBodyCapacity * 4L);
            }
        }

        internal void EnsureDampCapacity(int chunks, int members, int bodies)
        {
            var rt = Runtime!;
            if (chunks > ChunkCapacity)
            {
                rt.ReleaseBuffer(ref ChunkBody); rt.ReleaseBuffer(ref ChunkFrom); rt.ReleaseBuffer(ref ChunkTo); rt.ReleaseBuffer(ref ChunkPartial);
                ChunkCapacity = Math.Max(chunks, ChunkCapacity * 2);
                ChunkBody = rt.CreateBuffer(ChunkCapacity * 4L);
                ChunkFrom = rt.CreateBuffer(ChunkCapacity * 4L);
                ChunkTo = rt.CreateBuffer(ChunkCapacity * 4L);
                ChunkPartial = rt.CreateBuffer(ChunkCapacity * 64L);
            }

            if (members > BodyMemberCapacity)
            {
                rt.ReleaseBuffer(ref BodyMembers);
                BodyMemberCapacity = Math.Max(members, BodyMemberCapacity * 2);
                BodyMembers = rt.CreateBuffer(BodyMemberCapacity * 4L);
            }

            if (bodies > RigidCapacity)
            {
                rt.ReleaseBuffer(ref BodyStart); rt.ReleaseBuffer(ref BodyChunkStart); rt.ReleaseBuffer(ref BodyRigid);
                RigidCapacity = Math.Max(bodies, RigidCapacity * 2);
                BodyStart = rt.CreateBuffer((RigidCapacity + 1) * 4L);
                BodyChunkStart = rt.CreateBuffer((RigidCapacity + 1) * 4L);
                BodyRigid = rt.CreateBuffer(RigidCapacity * 36L);
            }
        }

        private readonly List<IntPtr> _particleBuffers = new();

        private IntPtr ParticleBuffer(long bytes)
        {
            var buffer = Runtime!.CreateBuffer(bytes);
            _particleBuffers.Add(buffer);
            return buffer;
        }

        private void ReleaseParticleBuffers()
        {
            for (int i = 0; i < _particleBuffers.Count; i++)
            {
                var buffer = _particleBuffers[i];
                Runtime!.ReleaseBuffer(ref buffer);
            }
            _particleBuffers.Clear();
            Capacity = 0;
        }
#endif

        /// <inheritdoc/>
        public void Dispose()
        {
#if !NETSTANDARD2_1
            Runtime?.Dispose();
            Runtime = null;
#endif
            _available = false;
        }
    }

    /// <summary>
    /// Compute backend types.
    /// </summary>
    public enum GpuBackendType
    {
        /// <summary>The multi-threaded CPU solver.</summary>
        CPU_SIMD,
        /// <summary>OpenCL (NVIDIA, AMD, Intel; deprecated but available on macOS).</summary>
        OpenCL,
        /// <summary>NVIDIA CUDA (not implemented: NVIDIA GPUs are used through OpenCL).</summary>
        CUDA,
        /// <summary>Apple Metal (not implemented: Macs are used through OpenCL).</summary>
        Metal,
        /// <summary>Vulkan compute (not implemented).</summary>
        Vulkan,
        /// <summary>DirectX 12 compute (not implemented).</summary>
        DirectCompute
    }
}
