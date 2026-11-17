using System;
using Artemis.Compute;

namespace Artemis.Particles
{
    // GPU path of the massive particle system.
    //
    // The particle state lives on the device while the simulation runs. It is copied
    // back only when the host actually looks at it (positions for rendering, a query,
    // an edit), and copied up again only after the host has changed it. The kernels
    // are a port of the CPU solver; the one feature that stays on the CPU is rigid
    // groups (shattering objects): while any exists, updates run on the CPU.
    public partial class MassiveParticleSystem
    {
        // Host and device copies of the particle state.
        private bool _hostStale;     // the device has newer data than the host arrays
        private bool _deviceStale = true; // the host has changed data the device lacks

        private int[] _intScratch = Array.Empty<int>();
        private readonly int[] _counterScratch = new int[4];

        /// <summary>
        /// Also accept OpenCL CPU devices (POCL, the Intel CPU runtime) when enabling
        /// <see cref="UseGPU"/>. Mainly useful for testing; set it before
        /// <see cref="UseGPU"/>.
        /// </summary>
        public bool AllowOpenCLCpuDevice { get; set; }

        /// <summary>
        /// Whether the last <see cref="Update"/> ran on the GPU.
        /// </summary>
        public bool LastUpdateUsedGpu { get; private set; }

        /// <summary>
        /// Why <see cref="UseGPU"/> could not be enabled, if it could not.
        /// </summary>
        public string? GpuError => _gpuCompute?.InitializationError ?? _gpuError;

        private string? _gpuError;

        /// <summary>
        /// Makes the host arrays current before they are read.
        /// </summary>
        private void SyncHost()
        {
            if (!_hostStale)
                return;
            _hostStale = false;
#if !NETSTANDARD2_1
            DownloadState();
#endif
        }

        /// <summary>
        /// Makes the host arrays current before they are changed, and remembers that
        /// the device copy has to be refreshed.
        /// </summary>
        private void TouchHost()
        {
            SyncHost();
            _deviceStale = true;
        }

        private bool CanUseGpu
#if NETSTANDARD2_1
            => false;
#else
            => _useGpu && _gpuCompute is { IsAvailable: true } && _aliveGroups == 0;
#endif

#if !NETSTANDARD2_1
        private const int LocalSize = 64;

        private void UpdateOnGpu(float dt, int subSteps)
        {
            var gpu = _gpuCompute!;
            var rt = gpu.Runtime!;
            int n = _count;

            gpu.EnsureCapacity(_capacity);
            if (_deviceStale)
            {
                UploadState();
                _deviceStale = false;
            }
            UploadSoftBodies();

            for (int step = 0; step < subSteps; step++)
                StepGpu(gpu, rt, n, dt);

            rt.Read(gpu.Counters, _counterScratch, 4);
            _activeCount = _counterScratch[1];
            _aliveCount = _counterScratch[2];
            _hostStale = true;
        }

        private void StepGpu(MassiveGpuCompute g, OpenCLRuntime rt, int n, float dt)
        {
            rt.FillInt(g.Counters, 0, 4);
            rt.Run("predict", n, LocalSize, n, dt, Gravity.X * dt, Gravity.Y * dt, Gravity.Z * dt, Damping,
                g.PX, g.PY, g.PZ, g.QX, g.QY, g.QZ, g.VX, g.VY, g.VZ,
                g.InvMass, g.Flags, g.Rest, g.Wake, g.DpX, g.DpY, g.DpZ, g.StX, g.StY, g.StZ, g.Counters);

            // Shape matching first, contacts after (see UpdateCpu).
            StepSoftBodiesGpu(g, rt, n, dt);

            rt.Read(g.Counters, _counterScratch, 1);
            _awakeColliders = _counterScratch[0];
            bool solved = EnableCollisions && _maxRadius > 0 && _awakeColliders > 0;

            if (solved)
            {
                float cell = Math.Max(2f * _maxRadius, 1e-4f);
                float invCell = 1f / cell;
                int table = 1024;
                while (table < n * 2 && table < (1 << 28))
                    table <<= 1;
                g.EnsureTable(table);
                int mask = table - 1;

                rt.FillInt(g.BucketCount, 0, table);
                rt.Run("assign_cells", n, LocalSize, n, invCell, mask, g.PX, g.PY, g.PZ, g.Flags, g.Cell, g.BucketCount);
                Scan(g, rt, g.BucketCount, g.BucketStart, table, 0);
                rt.Run("finalize_starts", 1, 1, g.BucketStart, g.BucketCount, table);
                rt.CopyInt(g.BucketStart, g.Cursor, table);
                rt.Run("scatter", n, LocalSize, n, g.Cell, g.Cursor, g.SortedIndex);
                rt.Run("sort_buckets", table, LocalSize, table, g.BucketStart, g.SortedIndex);
                rt.Run("gather", n, LocalSize, mask, invCell, g.BucketStart, g.SortedIndex,
                    g.PX, g.PY, g.PZ, g.QX, g.QY, g.QZ, g.Radius, g.InvMass, g.Flags, g.Body,
                    g.SX, g.SY, g.SZ, g.SR, g.SW, g.SMX, g.SMY, g.SMZ, g.SKey, g.SBody);

                GetContactParameters(dt, out float friction, out float relaxation, out float movingThresholdSq,
                    out float breakDistance, out bool canBreak, out float stackingBias, out float upX, out float upY, out float upZ);

                rt.Run("solve", n, LocalSize, mask, invCell, friction, relaxation, movingThresholdSq,
                    canBreak ? breakDistance : 0f, canBreak ? 1 : 0, stackingBias, upX, upY, upZ,
                    g.BucketStart, g.SortedIndex, g.SX, g.SY, g.SZ, g.SR, g.SW, g.SMX, g.SMY, g.SMZ, g.SKey, g.SBody,
                    g.Wake, g.DpX, g.DpY, g.DpZ, g.StX, g.StY, g.StZ);
                rt.Run("apply_contacts", n, LocalSize, mask, g.BucketStart, g.SortedIndex, g.SW, g.PX, g.PY, g.PZ, g.DpX, g.DpY, g.DpZ);
            }

            rt.Run("finish", n, LocalSize, n, dt,
                (float)Bounds.Min.X, (float)Bounds.Min.Y, (float)Bounds.Min.Z,
                (float)Bounds.Max.X, (float)Bounds.Max.Y, (float)Bounds.Max.Z,
                Restitution, Math.Max(Friction, 0), AllowSleeping ? 1 : 0, SleepSpeed * SleepSpeed,
                Math.Clamp(SleepSubSteps, 1, 255), solved ? 1 : 0,
                g.PX, g.PY, g.PZ, g.QX, g.QY, g.QZ, g.VX, g.VY, g.VZ,
                g.InvMass, g.Radius, g.Lifetime, g.Flags, g.Rest, g.Body,
                g.DpX, g.DpY, g.DpZ, g.StX, g.StY, g.StZ, g.Counters);

            DampSoftBodiesGpu(g, rt, n);
        }

        /// <summary>
        /// Exclusive prefix sum (recursive over blocks of 1024).
        /// </summary>
        private static void Scan(MassiveGpuCompute g, OpenCLRuntime rt, IntPtr input, IntPtr output, int length, int level)
        {
            int blocks = (length + ParticleKernels.ScanBlock - 1) / ParticleKernels.ScanBlock;
            var sums = g.ScanSums[level];
            rt.Run("scan_block", blocks * ParticleKernels.ScanGroupSize, ParticleKernels.ScanGroupSize, input, output, sums, length);
            if (blocks > 1)
            {
                Scan(g, rt, sums, sums, blocks, level + 1);
                rt.Run("add_offsets", length, LocalSize, output, sums, length);
            }
        }

        private int[] IntScratch(int length)
        {
            if (_intScratch.Length < length)
                _intScratch = new int[Math.Max(length, 1024)];
            return _intScratch;
        }

        private void UploadState()
        {
            var g = _gpuCompute!;
            var rt = g.Runtime!;
            int n = _count;

            rt.Write(g.PX, _posX, n); rt.Write(g.PY, _posY, n); rt.Write(g.PZ, _posZ, n);
            rt.Write(g.VX, _velX, n); rt.Write(g.VY, _velY, n); rt.Write(g.VZ, _velZ, n);
            rt.Write(g.InvMass, _invMass, n);
            rt.Write(g.Radius, _radius, n);
            rt.Write(g.Lifetime, _lifetime, n);
            rt.Write(g.Flags, _flags, n);
            rt.Write(g.Rest, _rest, n);

            var scratch = IntScratch(n);
            for (int i = 0; i < n; i++) scratch[i] = _wake[i];
            rt.Write(g.Wake, scratch, n);

            for (int i = 0; i < n; i++) scratch[i] = BodyTag(i);
            rt.Write(g.Body, scratch, n);
        }

        private void DownloadState()
        {
            var g = _gpuCompute;
            if (g?.Runtime == null)
                return;
            var rt = g.Runtime;
            int n = _count;

            rt.Read(g.PX, _posX, n); rt.Read(g.PY, _posY, n); rt.Read(g.PZ, _posZ, n);
            rt.Read(g.VX, _velX, n); rt.Read(g.VY, _velY, n); rt.Read(g.VZ, _velZ, n);
            rt.Read(g.Lifetime, _lifetime, n);
            rt.Read(g.Flags, _flags, n);
            rt.Read(g.Rest, _rest, n);

            var scratch = IntScratch(n);
            rt.Read(g.Wake, scratch, n);
            for (int i = 0; i < n; i++) _wake[i] = (byte)scratch[i];

            DownloadSoftBodies();
        }
#endif

        /// <summary>
        /// Id used to keep particles of the same object from colliding with each
        /// other: the rigid group, the tagged soft body, or -1 for loose particles.
        /// </summary>
        private int BodyTag(int i)
        {
            if (_groupOf[i] >= 0) return _groupOf[i];
            if (_softOf[i] >= 0) return _softOf[i] | SoftBodyTag;
            return -1;
        }

        private const int SoftBodyTag = 0x40000000;

        /// <summary>
        /// Contact parameters shared by the CPU and GPU solvers.
        /// </summary>
        private void GetContactParameters(float dt, out float friction, out float relaxation, out float movingThresholdSq,
            out float breakDistance, out bool canBreak, out float stackingBias, out float upX, out float upY, out float upZ)
        {
            friction = Math.Max(Friction, 0);
            relaxation = ContactRelaxation;
            movingThresholdSq = WakeSpeed * dt * WakeSpeed * dt;
            breakDistance = FrozenBreakSpeed * dt;
            canBreak = !float.IsInfinity(FrozenBreakSpeed);
            stackingBias = Math.Clamp(StackingBias, 0f, 1f);
            upX = upY = upZ = 0;

            float gravityLength = Gravity.Length();
            if (gravityLength > 1e-6f)
            {
                upX = -Gravity.X / gravityLength;
                upY = -Gravity.Y / gravityLength;
                upZ = -Gravity.Z / gravityLength;
            }
            else
            {
                stackingBias = 0;
            }
        }
    }
}
