using System;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Runtime.Intrinsics.Arm;
using System.Threading.Tasks;

namespace Artemis.Compute
{
    /// <summary>
    /// High-performance SIMD helper methods for physics calculations.
    /// Supports SSE, AVX, AVX-512 (x64) and NEON (ARM64/Apple Silicon).
    /// Auto-detects the best SIMD width for the current processor.
    /// </summary>
    public static class SimdHelpers
    {
        #region SIMD Capability Detection

        /// <summary>
        /// Whether AVX-512 is supported.
        /// </summary>
        public static readonly bool HasAvx512 = Avx512F.IsSupported;

        /// <summary>
        /// Whether AVX2 is supported.
        /// </summary>
        public static readonly bool HasAvx2 = Avx2.IsSupported;

        /// <summary>
        /// Whether AVX is supported.
        /// </summary>
        public static readonly bool HasAvx = Avx.IsSupported;

        /// <summary>
        /// Whether SSE is supported.
        /// </summary>
        public static readonly bool HasSse = Sse.IsSupported;

        /// <summary>
        /// Whether ARM NEON is supported (Apple Silicon M1/M2/M3).
        /// </summary>
        public static readonly bool HasNeon = AdvSimd.IsSupported;

        /// <summary>
        /// Whether ARM64 advanced SIMD is supported.
        /// </summary>
        public static readonly bool HasArm64AdvSimd = AdvSimd.Arm64.IsSupported;

        /// <summary>
        /// Optimal SIMD width in floats for this processor.
        /// </summary>
        public static readonly int OptimalSimdWidth =
            HasAvx512 ? 16 :
            HasAvx2 || HasAvx ? 8 :
            HasSse || HasNeon ? 4 :
            Vector<float>.Count;

        /// <summary>
        /// Gets a description of the available SIMD features.
        /// </summary>
        public static string SimdCapabilities
        {
            get
            {
                if (HasAvx512) return "AVX-512 (512-bit, 16 floats)";
                if (HasAvx2) return "AVX2 (256-bit, 8 floats)";
                if (HasAvx) return "AVX (256-bit, 8 floats)";
                if (HasNeon && HasArm64AdvSimd) return "ARM64 NEON (128-bit, 4 floats) - Apple Silicon";
                if (HasNeon) return "ARM NEON (128-bit, 4 floats)";
                if (HasSse) return "SSE (128-bit, 4 floats)";
                return $"Generic SIMD ({Vector<float>.Count} floats)";
            }
        }

        #endregion

        #region Vector Addition (Float Arrays)

        /// <summary>
        /// Adds two float arrays element-wise using the best available SIMD.
        /// result[i] = a[i] + b[i]
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void Add(float[] a, float[] b, float[] result, int count)
        {
            int i = 0;

            if (HasAvx512 && count >= 16)
            {
                AddAvx512(a, b, result, count, ref i);
            }
            else if (HasAvx && count >= 8)
            {
                AddAvx(a, b, result, count, ref i);
            }
            else if (HasNeon && count >= 4)
            {
                AddNeon(a, b, result, count, ref i);
            }
            else if (HasSse && count >= 4)
            {
                AddSse(a, b, result, count, ref i);
            }

            // Scalar remainder
            for (; i < count; i++)
            {
                result[i] = a[i] + b[i];
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AddAvx512(float[] a, float[] b, float[] result, int count, ref int i)
        {
            for (; i <= count - 16; i += 16)
            {
                var va = Avx512F.LoadVector512(ref a[i]);
                var vb = Avx512F.LoadVector512(ref b[i]);
                var vr = Avx512F.Add(va, vb);
                Avx512F.Store(ref result[i], vr);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AddAvx(float[] a, float[] b, float[] result, int count, ref int i)
        {
            for (; i <= count - 8; i += 8)
            {
                var va = Avx.LoadVector256(ref a[i]);
                var vb = Avx.LoadVector256(ref b[i]);
                var vr = Avx.Add(va, vb);
                Avx.Store(ref result[i], vr);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AddSse(float[] a, float[] b, float[] result, int count, ref int i)
        {
            for (; i <= count - 4; i += 4)
            {
                var va = Sse.LoadVector128(ref a[i]);
                var vb = Sse.LoadVector128(ref b[i]);
                var vr = Sse.Add(va, vb);
                Sse.Store(ref result[i], vr);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void AddNeon(float[] a, float[] b, float[] result, int count, ref int i)
        {
            for (; i <= count - 4; i += 4)
            {
                var va = AdvSimd.LoadVector128(ref a[i]);
                var vb = AdvSimd.LoadVector128(ref b[i]);
                var vr = AdvSimd.Add(va, vb);
                AdvSimd.Store(ref result[i], vr);
            }
        }

        #endregion

        #region Multiply-Add (FMA) Operations

        /// <summary>
        /// Performs fused multiply-add: result[i] = a[i] * b + c[i]
        /// Uses FMA instructions when available for better precision and performance.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void MultiplyAdd(float[] a, float scalar, float[] c, float[] result, int count)
        {
            int i = 0;

            if (Fma.IsSupported && count >= 8)
            {
                var scalarVec = Vector256.Create(scalar);
                for (; i <= count - 8; i += 8)
                {
                    var va = Avx.LoadVector256(ref a[i]);
                    var vc = Avx.LoadVector256(ref c[i]);
                    var vr = Fma.MultiplyAdd(va, scalarVec, vc);
                    Avx.Store(ref result[i], vr);
                }
            }
            else if (HasNeon && HasArm64AdvSimd && count >= 4)
            {
                var scalarVec = Vector128.Create(scalar);
                for (; i <= count - 4; i += 4)
                {
                    var va = AdvSimd.LoadVector128(ref a[i]);
                    var vc = AdvSimd.LoadVector128(ref c[i]);
                    // ARM64 has FMA
                    var vr = AdvSimd.Arm64.FusedMultiplyAdd(vc, va, scalarVec);
                    AdvSimd.Store(ref result[i], vr);
                }
            }

            // Scalar remainder
            for (; i < count; i++)
            {
                result[i] = a[i] * scalar + c[i];
            }
        }

        /// <summary>
        /// In-place multiply-add: arr[i] = arr[i] * scalar + addend
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void MultiplyAddInPlace(float[] arr, float scalar, float addend, int count)
        {
            int i = 0;

            if (Fma.IsSupported && count >= 8)
            {
                var scalarVec = Vector256.Create(scalar);
                var addVec = Vector256.Create(addend);
                for (; i <= count - 8; i += 8)
                {
                    var va = Avx.LoadVector256(ref arr[i]);
                    var vr = Fma.MultiplyAdd(va, scalarVec, addVec);
                    Avx.Store(ref arr[i], vr);
                }
            }
            else if (HasNeon && count >= 4)
            {
                var scalarVec = Vector128.Create(scalar);
                var addVec = Vector128.Create(addend);
                for (; i <= count - 4; i += 4)
                {
                    var va = AdvSimd.LoadVector128(ref arr[i]);
                    var vr = AdvSimd.Add(AdvSimd.Multiply(va, scalarVec), addVec);
                    AdvSimd.Store(ref arr[i], vr);
                }
            }

            for (; i < count; i++)
            {
                arr[i] = arr[i] * scalar + addend;
            }
        }

        #endregion

        #region Distance Calculations

        /// <summary>
        /// Computes squared distances between particle pairs.
        /// Optimized for spatial hash collision detection.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void ComputeDistancesSq(
            float[] x1, float[] y1, float[] z1,
            float[] x2, float[] y2, float[] z2,
            float[] distSq, int count)
        {
            int i = 0;

            if (HasAvx && count >= 8)
            {
                for (; i <= count - 8; i += 8)
                {
                    var dx = Avx.Subtract(Avx.LoadVector256(ref x2[i]), Avx.LoadVector256(ref x1[i]));
                    var dy = Avx.Subtract(Avx.LoadVector256(ref y2[i]), Avx.LoadVector256(ref y1[i]));
                    var dz = Avx.Subtract(Avx.LoadVector256(ref z2[i]), Avx.LoadVector256(ref z1[i]));

                    var dx2 = Avx.Multiply(dx, dx);
                    var dy2 = Avx.Multiply(dy, dy);
                    var dz2 = Avx.Multiply(dz, dz);

                    var result = Avx.Add(Avx.Add(dx2, dy2), dz2);
                    Avx.Store(ref distSq[i], result);
                }
            }
            else if (HasNeon && count >= 4)
            {
                for (; i <= count - 4; i += 4)
                {
                    var dx = AdvSimd.Subtract(AdvSimd.LoadVector128(ref x2[i]), AdvSimd.LoadVector128(ref x1[i]));
                    var dy = AdvSimd.Subtract(AdvSimd.LoadVector128(ref y2[i]), AdvSimd.LoadVector128(ref y1[i]));
                    var dz = AdvSimd.Subtract(AdvSimd.LoadVector128(ref z2[i]), AdvSimd.LoadVector128(ref z1[i]));

                    var dx2 = AdvSimd.Multiply(dx, dx);
                    var dy2 = AdvSimd.Multiply(dy, dy);
                    var dz2 = AdvSimd.Multiply(dz, dz);

                    var result = AdvSimd.Add(AdvSimd.Add(dx2, dy2), dz2);
                    AdvSimd.Store(ref distSq[i], result);
                }
            }

            for (; i < count; i++)
            {
                float dx = x2[i] - x1[i];
                float dy = y2[i] - y1[i];
                float dz = z2[i] - z1[i];
                distSq[i] = dx * dx + dy * dy + dz * dz;
            }
        }

        #endregion

        #region Particle Integration

        /// <summary>
        /// High-performance particle velocity integration.
        /// Updates velocity arrays with gravity and damping.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void IntegrateVelocities(
            float[] velX, float[] velY, float[] velZ,
            float gravityX, float gravityY, float gravityZ,
            float damping, float dt, int count)
        {
            float gxDt = gravityX * dt;
            float gyDt = gravityY * dt;
            float gzDt = gravityZ * dt;

            int i = 0;

            if (HasAvx && count >= 8)
            {
                var gxVec = Vector256.Create(gxDt);
                var gyVec = Vector256.Create(gyDt);
                var gzVec = Vector256.Create(gzDt);
                var dampVec = Vector256.Create(damping);

                for (; i <= count - 8; i += 8)
                {
                    var vx = Avx.LoadVector256(ref velX[i]);
                    var vy = Avx.LoadVector256(ref velY[i]);
                    var vz = Avx.LoadVector256(ref velZ[i]);

                    vx = Avx.Add(vx, gxVec);
                    vy = Avx.Add(vy, gyVec);
                    vz = Avx.Add(vz, gzVec);

                    vx = Avx.Multiply(vx, dampVec);
                    vy = Avx.Multiply(vy, dampVec);
                    vz = Avx.Multiply(vz, dampVec);

                    Avx.Store(ref velX[i], vx);
                    Avx.Store(ref velY[i], vy);
                    Avx.Store(ref velZ[i], vz);
                }
            }
            else if (HasNeon && count >= 4)
            {
                var gxVec = Vector128.Create(gxDt);
                var gyVec = Vector128.Create(gyDt);
                var gzVec = Vector128.Create(gzDt);
                var dampVec = Vector128.Create(damping);

                for (; i <= count - 4; i += 4)
                {
                    var vx = AdvSimd.LoadVector128(ref velX[i]);
                    var vy = AdvSimd.LoadVector128(ref velY[i]);
                    var vz = AdvSimd.LoadVector128(ref velZ[i]);

                    vx = AdvSimd.Add(vx, gxVec);
                    vy = AdvSimd.Add(vy, gyVec);
                    vz = AdvSimd.Add(vz, gzVec);

                    vx = AdvSimd.Multiply(vx, dampVec);
                    vy = AdvSimd.Multiply(vy, dampVec);
                    vz = AdvSimd.Multiply(vz, dampVec);

                    AdvSimd.Store(ref velX[i], vx);
                    AdvSimd.Store(ref velY[i], vy);
                    AdvSimd.Store(ref velZ[i], vz);
                }
            }

            for (; i < count; i++)
            {
                velX[i] = (velX[i] + gxDt) * damping;
                velY[i] = (velY[i] + gyDt) * damping;
                velZ[i] = (velZ[i] + gzDt) * damping;
            }
        }

        /// <summary>
        /// High-performance particle position integration.
        /// Updates position arrays based on velocity.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void IntegratePositions(
            float[] posX, float[] posY, float[] posZ,
            float[] velX, float[] velY, float[] velZ,
            float dt, int count)
        {
            int i = 0;

            if (Fma.IsSupported && count >= 8)
            {
                var dtVec = Vector256.Create(dt);

                for (; i <= count - 8; i += 8)
                {
                    var px = Avx.LoadVector256(ref posX[i]);
                    var py = Avx.LoadVector256(ref posY[i]);
                    var pz = Avx.LoadVector256(ref posZ[i]);
                    var vx = Avx.LoadVector256(ref velX[i]);
                    var vy = Avx.LoadVector256(ref velY[i]);
                    var vz = Avx.LoadVector256(ref velZ[i]);

                    px = Fma.MultiplyAdd(vx, dtVec, px);
                    py = Fma.MultiplyAdd(vy, dtVec, py);
                    pz = Fma.MultiplyAdd(vz, dtVec, pz);

                    Avx.Store(ref posX[i], px);
                    Avx.Store(ref posY[i], py);
                    Avx.Store(ref posZ[i], pz);
                }
            }
            else if (HasNeon && HasArm64AdvSimd && count >= 4)
            {
                var dtVec = Vector128.Create(dt);

                for (; i <= count - 4; i += 4)
                {
                    var px = AdvSimd.LoadVector128(ref posX[i]);
                    var py = AdvSimd.LoadVector128(ref posY[i]);
                    var pz = AdvSimd.LoadVector128(ref posZ[i]);
                    var vx = AdvSimd.LoadVector128(ref velX[i]);
                    var vy = AdvSimd.LoadVector128(ref velY[i]);
                    var vz = AdvSimd.LoadVector128(ref velZ[i]);

                    // ARM64 FMA: result = c + a * b
                    px = AdvSimd.Arm64.FusedMultiplyAdd(px, vx, dtVec);
                    py = AdvSimd.Arm64.FusedMultiplyAdd(py, vy, dtVec);
                    pz = AdvSimd.Arm64.FusedMultiplyAdd(pz, vz, dtVec);

                    AdvSimd.Store(ref posX[i], px);
                    AdvSimd.Store(ref posY[i], py);
                    AdvSimd.Store(ref posZ[i], pz);
                }
            }

            for (; i < count; i++)
            {
                posX[i] += velX[i] * dt;
                posY[i] += velY[i] * dt;
                posZ[i] += velZ[i] * dt;
            }
        }

        #endregion

        #region SPH Kernel Calculations

        /// <summary>
        /// Computes the Poly6 kernel for SPH density calculation.
        /// W(r, h) = 315 / (64 * pi * h^9) * (h^2 - r^2)^3
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void ComputePoly6Kernel(
            float[] distSq, float[] result, float hSq, float poly6Coeff, int count)
        {
            int i = 0;

            if (HasAvx && count >= 8)
            {
                var hSqVec = Vector256.Create(hSq);
                var coeffVec = Vector256.Create(poly6Coeff);
                var zeroVec = Vector256<float>.Zero;

                for (; i <= count - 8; i += 8)
                {
                    var r2 = Avx.LoadVector256(ref distSq[i]);
                    var diff = Avx.Subtract(hSqVec, r2);

                    // Only compute if within radius (diff > 0)
                    var mask = Avx.Compare(diff, zeroVec, FloatComparisonMode.OrderedGreaterThanSignaling);
                    diff = Avx.And(diff, mask);

                    // (h^2 - r^2)^3
                    var diff2 = Avx.Multiply(diff, diff);
                    var diff3 = Avx.Multiply(diff2, diff);

                    var w = Avx.Multiply(diff3, coeffVec);
                    Avx.Store(ref result[i], w);
                }
            }
            else if (HasNeon && count >= 4)
            {
                var hSqVec = Vector128.Create(hSq);
                var coeffVec = Vector128.Create(poly6Coeff);

                for (; i <= count - 4; i += 4)
                {
                    var r2 = AdvSimd.LoadVector128(ref distSq[i]);
                    var diff = AdvSimd.Subtract(hSqVec, r2);

                    // Clamp negative to zero
                    diff = AdvSimd.Max(diff, Vector128<float>.Zero);

                    var diff2 = AdvSimd.Multiply(diff, diff);
                    var diff3 = AdvSimd.Multiply(diff2, diff);

                    var w = AdvSimd.Multiply(diff3, coeffVec);
                    AdvSimd.Store(ref result[i], w);
                }
            }

            for (; i < count; i++)
            {
                float diff = hSq - distSq[i];
                if (diff > 0)
                {
                    result[i] = poly6Coeff * diff * diff * diff;
                }
                else
                {
                    result[i] = 0;
                }
            }
        }

        #endregion

        #region Parallel SIMD Operations

        /// <summary>
        /// Parallel SIMD integration for massive particle counts.
        /// Splits work across multiple threads, each using SIMD.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public static void ParallelIntegrate(
            float[] posX, float[] posY, float[] posZ,
            float[] velX, float[] velY, float[] velZ,
            float gravityX, float gravityY, float gravityZ,
            float damping, float dt, int count,
            int threadCount = 0)
        {
            if (threadCount <= 0)
                threadCount = Environment.ProcessorCount;

            int chunkSize = (count + threadCount - 1) / threadCount;

            Parallel.For(0, threadCount, new ParallelOptions { MaxDegreeOfParallelism = threadCount }, threadId =>
            {
                int start = threadId * chunkSize;
                int end = Math.Min(start + chunkSize, count);
                int localCount = end - start;

                if (localCount <= 0) return;

                // Create spans for this chunk
                var pxSpan = posX.AsSpan(start, localCount);
                var pySpan = posY.AsSpan(start, localCount);
                var pzSpan = posZ.AsSpan(start, localCount);
                var vxSpan = velX.AsSpan(start, localCount);
                var vySpan = velY.AsSpan(start, localCount);
                var vzSpan = velZ.AsSpan(start, localCount);

                // Integrate this chunk
                IntegrateChunk(
                    ref MemoryMarshal.GetReference(pxSpan),
                    ref MemoryMarshal.GetReference(pySpan),
                    ref MemoryMarshal.GetReference(pzSpan),
                    ref MemoryMarshal.GetReference(vxSpan),
                    ref MemoryMarshal.GetReference(vySpan),
                    ref MemoryMarshal.GetReference(vzSpan),
                    gravityX, gravityY, gravityZ,
                    damping, dt, localCount
                );
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        private static void IntegrateChunk(
            ref float posX, ref float posY, ref float posZ,
            ref float velX, ref float velY, ref float velZ,
            float gravityX, float gravityY, float gravityZ,
            float damping, float dt, int count)
        {
            float gxDt = gravityX * dt;
            float gyDt = gravityY * dt;
            float gzDt = gravityZ * dt;

            int i = 0;

            if (HasAvx && count >= 8)
            {
                var gxVec = Vector256.Create(gxDt);
                var gyVec = Vector256.Create(gyDt);
                var gzVec = Vector256.Create(gzDt);
                var dampVec = Vector256.Create(damping);
                var dtVec = Vector256.Create(dt);

                for (; i <= count - 8; i += 8)
                {
                    // Load velocities
                    var vx = Avx.LoadVector256(ref Unsafe.Add(ref velX, i));
                    var vy = Avx.LoadVector256(ref Unsafe.Add(ref velY, i));
                    var vz = Avx.LoadVector256(ref Unsafe.Add(ref velZ, i));

                    // Apply gravity
                    vx = Avx.Add(vx, gxVec);
                    vy = Avx.Add(vy, gyVec);
                    vz = Avx.Add(vz, gzVec);

                    // Apply damping
                    vx = Avx.Multiply(vx, dampVec);
                    vy = Avx.Multiply(vy, dampVec);
                    vz = Avx.Multiply(vz, dampVec);

                    // Store velocities
                    Avx.Store(ref Unsafe.Add(ref velX, i), vx);
                    Avx.Store(ref Unsafe.Add(ref velY, i), vy);
                    Avx.Store(ref Unsafe.Add(ref velZ, i), vz);

                    // Load positions
                    var px = Avx.LoadVector256(ref Unsafe.Add(ref posX, i));
                    var py = Avx.LoadVector256(ref Unsafe.Add(ref posY, i));
                    var pz = Avx.LoadVector256(ref Unsafe.Add(ref posZ, i));

                    // Update positions: p += v * dt
                    if (Fma.IsSupported)
                    {
                        px = Fma.MultiplyAdd(vx, dtVec, px);
                        py = Fma.MultiplyAdd(vy, dtVec, py);
                        pz = Fma.MultiplyAdd(vz, dtVec, pz);
                    }
                    else
                    {
                        px = Avx.Add(px, Avx.Multiply(vx, dtVec));
                        py = Avx.Add(py, Avx.Multiply(vy, dtVec));
                        pz = Avx.Add(pz, Avx.Multiply(vz, dtVec));
                    }

                    // Store positions
                    Avx.Store(ref Unsafe.Add(ref posX, i), px);
                    Avx.Store(ref Unsafe.Add(ref posY, i), py);
                    Avx.Store(ref Unsafe.Add(ref posZ, i), pz);
                }
            }
            else if (HasNeon && count >= 4)
            {
                var gxVec = Vector128.Create(gxDt);
                var gyVec = Vector128.Create(gyDt);
                var gzVec = Vector128.Create(gzDt);
                var dampVec = Vector128.Create(damping);
                var dtVec = Vector128.Create(dt);

                for (; i <= count - 4; i += 4)
                {
                    var vx = AdvSimd.LoadVector128(ref Unsafe.Add(ref velX, i));
                    var vy = AdvSimd.LoadVector128(ref Unsafe.Add(ref velY, i));
                    var vz = AdvSimd.LoadVector128(ref Unsafe.Add(ref velZ, i));

                    vx = AdvSimd.Add(vx, gxVec);
                    vy = AdvSimd.Add(vy, gyVec);
                    vz = AdvSimd.Add(vz, gzVec);

                    vx = AdvSimd.Multiply(vx, dampVec);
                    vy = AdvSimd.Multiply(vy, dampVec);
                    vz = AdvSimd.Multiply(vz, dampVec);

                    AdvSimd.Store(ref Unsafe.Add(ref velX, i), vx);
                    AdvSimd.Store(ref Unsafe.Add(ref velY, i), vy);
                    AdvSimd.Store(ref Unsafe.Add(ref velZ, i), vz);

                    var px = AdvSimd.LoadVector128(ref Unsafe.Add(ref posX, i));
                    var py = AdvSimd.LoadVector128(ref Unsafe.Add(ref posY, i));
                    var pz = AdvSimd.LoadVector128(ref Unsafe.Add(ref posZ, i));

                    if (HasArm64AdvSimd)
                    {
                        px = AdvSimd.Arm64.FusedMultiplyAdd(px, vx, dtVec);
                        py = AdvSimd.Arm64.FusedMultiplyAdd(py, vy, dtVec);
                        pz = AdvSimd.Arm64.FusedMultiplyAdd(pz, vz, dtVec);
                    }
                    else
                    {
                        px = AdvSimd.Add(px, AdvSimd.Multiply(vx, dtVec));
                        py = AdvSimd.Add(py, AdvSimd.Multiply(vy, dtVec));
                        pz = AdvSimd.Add(pz, AdvSimd.Multiply(vz, dtVec));
                    }

                    AdvSimd.Store(ref Unsafe.Add(ref posX, i), px);
                    AdvSimd.Store(ref Unsafe.Add(ref posY, i), py);
                    AdvSimd.Store(ref Unsafe.Add(ref posZ, i), pz);
                }
            }

            // Scalar remainder
            for (; i < count; i++)
            {
                ref float vx = ref Unsafe.Add(ref velX, i);
                ref float vy = ref Unsafe.Add(ref velY, i);
                ref float vz = ref Unsafe.Add(ref velZ, i);

                vx = (vx + gxDt) * damping;
                vy = (vy + gyDt) * damping;
                vz = (vz + gzDt) * damping;

                Unsafe.Add(ref posX, i) += vx * dt;
                Unsafe.Add(ref posY, i) += vy * dt;
                Unsafe.Add(ref posZ, i) += vz * dt;
            }
        }

        #endregion
    }
}
