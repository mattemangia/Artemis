#if !NETSTANDARD2_1
using System;
using System.Collections.Generic;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;

namespace Artemis.Compute
{
    /// <summary>
    /// Thrown when an OpenCL call fails.
    /// </summary>
    public sealed class OpenCLException : Exception
    {
        /// <summary>OpenCL error code (negative).</summary>
        public int ErrorCode { get; }

        internal OpenCLException(string call, int code, string? details = null)
            : base($"{call} failed: {OpenCL.ErrorName(code)} ({code})" + (details != null ? Environment.NewLine + details : ""))
        {
            ErrorCode = code;
        }
    }

    /// <summary>
    /// Description of an OpenCL compute device.
    /// </summary>
    public sealed class OpenCLDeviceInfo
    {
        internal IntPtr Platform;
        internal IntPtr Device;

        /// <summary>Device name as reported by the driver.</summary>
        public string Name { get; internal set; } = "";

        /// <summary>Device vendor.</summary>
        public string Vendor { get; internal set; } = "";

        /// <summary>Name of the OpenCL platform (driver) the device belongs to.</summary>
        public string PlatformName { get; internal set; } = "";

        /// <summary>OpenCL version string of the device.</summary>
        public string Version { get; internal set; } = "";

        /// <summary>Whether this is a GPU.</summary>
        public bool IsGpu { get; internal set; }

        /// <summary>Whether this is a CPU device (for example POCL or the Intel CPU runtime).</summary>
        public bool IsCpu { get; internal set; }

        /// <summary>Number of compute units.</summary>
        public int ComputeUnits { get; internal set; }

        /// <summary>Global memory in bytes.</summary>
        public long GlobalMemory { get; internal set; }

        /// <summary>Largest single allocation in bytes.</summary>
        public long MaxAllocation { get; internal set; }

        /// <summary>Largest work-group size.</summary>
        public int MaxWorkGroupSize { get; internal set; }

        /// <inheritdoc/>
        public override string ToString()
            => $"{Name} ({(IsGpu ? "GPU" : IsCpu ? "CPU" : "accelerator")}, {Vendor}, {ComputeUnits} CUs, {GlobalMemory / (1024 * 1024)} MB, {PlatformName})";
    }

    /// <summary>
    /// Minimal OpenCL 1.2 bindings, loaded at run time. Nothing here requires an
    /// OpenCL SDK at build time: on machines without a driver
    /// <see cref="IsAvailable"/> is simply false.
    /// </summary>
    public static class OpenCL
    {
        private const string Library = "OpenCL";

        #region Constants

        internal const int Success = 0;
        internal const ulong DeviceTypeCpu = 1 << 1;
        internal const ulong DeviceTypeGpu = 1 << 2;
        internal const ulong DeviceTypeAccelerator = 1 << 3;
        internal const ulong DeviceTypeAll = 0xFFFFFFFF;

        private const uint PlatformName = 0x0902;
        private const uint DeviceType = 0x1000;
        private const uint DeviceMaxComputeUnits = 0x1002;
        private const uint DeviceMaxWorkGroupSize = 0x1004;
        private const uint DeviceMaxMemAllocSize = 0x1010;
        private const uint DeviceGlobalMemSize = 0x101F;
        private const uint DeviceAvailable = 0x1027;
        private const uint DeviceName = 0x102B;
        private const uint DeviceVendor = 0x102C;
        private const uint DeviceVersion = 0x102F;
        internal const uint ProgramBuildLog = 0x1183;
        internal const ulong MemReadWrite = 1;

        #endregion

        #region Loading

        private static readonly object Gate = new();
        private static bool _probed;
        private static bool _available;
        private static List<OpenCLDeviceInfo>? _devices;

        static OpenCL()
        {
            try
            {
                NativeLibrary.SetDllImportResolver(typeof(OpenCL).Assembly, Resolve);
            }
            catch (InvalidOperationException)
            {
                // A resolver is already registered for this assembly.
            }
        }

        private static IntPtr Resolve(string name, Assembly assembly, DllImportSearchPath? path)
        {
            if (name != Library)
                return IntPtr.Zero;

            string[] candidates =
                RuntimeInformation.IsOSPlatform(OSPlatform.Windows) ? new[] { "OpenCL.dll" } :
                RuntimeInformation.IsOSPlatform(OSPlatform.OSX) ? new[] { "/System/Library/Frameworks/OpenCL.framework/OpenCL" } :
                new[] { "libOpenCL.so.1", "libOpenCL.so" };

            foreach (var candidate in candidates)
            {
                if (NativeLibrary.TryLoad(candidate, assembly, path, out var handle))
                    return handle;
            }

            return IntPtr.Zero;
        }

        /// <summary>
        /// Whether an OpenCL driver with at least one usable device is installed.
        /// </summary>
        public static bool IsAvailable
        {
            get
            {
                Probe();
                return _available;
            }
        }

        /// <summary>
        /// All OpenCL devices of all installed platforms.
        /// </summary>
        public static IReadOnlyList<OpenCLDeviceInfo> Devices
        {
            get
            {
                Probe();
                return _devices!;
            }
        }

        /// <summary>
        /// The device to use: the GPU with the most compute units, then other
        /// accelerators, then (only when <paramref name="allowCpu"/> is set) CPU devices.
        /// </summary>
        public static OpenCLDeviceInfo? SelectDevice(bool allowCpu = false)
        {
            OpenCLDeviceInfo? best = null;
            int Rank(OpenCLDeviceInfo d) => d.IsGpu ? 2 : d.IsCpu ? 0 : 1;

            foreach (var device in Devices)
            {
                if (device.IsCpu && !allowCpu)
                    continue;

                if (best == null || Rank(device) > Rank(best) ||
                    (Rank(device) == Rank(best) && device.ComputeUnits > best.ComputeUnits))
                    best = device;
            }

            return best;
        }

        private static void Probe()
        {
            lock (Gate)
            {
                if (_probed)
                    return;
                _probed = true;
                _devices = new List<OpenCLDeviceInfo>();

                try
                {
                    if (clGetPlatformIDs(0, null, out uint platformCount) != Success || platformCount == 0)
                        return;

                    var platforms = new IntPtr[platformCount];
                    if (clGetPlatformIDs(platformCount, platforms, out _) != Success)
                        return;

                    foreach (var platform in platforms)
                    {
                        if (clGetDeviceIDs(platform, DeviceTypeAll, 0, null, out uint deviceCount) != Success || deviceCount == 0)
                            continue;

                        var devices = new IntPtr[deviceCount];
                        if (clGetDeviceIDs(platform, DeviceTypeAll, deviceCount, devices, out _) != Success)
                            continue;

                        string platformName = GetPlatformString(platform, PlatformName);
                        foreach (var device in devices)
                        {
                            if (GetDeviceUInt(device, DeviceAvailable) == 0)
                                continue;

                            ulong type = GetDeviceULong(device, DeviceType);
                            _devices.Add(new OpenCLDeviceInfo
                            {
                                Platform = platform,
                                Device = device,
                                Name = GetDeviceString(device, DeviceName).Trim(),
                                Vendor = GetDeviceString(device, DeviceVendor).Trim(),
                                Version = GetDeviceString(device, DeviceVersion).Trim(),
                                PlatformName = platformName.Trim(),
                                IsGpu = (type & DeviceTypeGpu) != 0,
                                IsCpu = (type & DeviceTypeCpu) != 0,
                                ComputeUnits = (int)GetDeviceUInt(device, DeviceMaxComputeUnits),
                                GlobalMemory = (long)GetDeviceULong(device, DeviceGlobalMemSize),
                                MaxAllocation = (long)GetDeviceULong(device, DeviceMaxMemAllocSize),
                                MaxWorkGroupSize = (int)GetDeviceSize(device, DeviceMaxWorkGroupSize)
                            });
                        }
                    }

                    _available = _devices.Count > 0;
                }
                catch (DllNotFoundException)
                {
                }
                catch (EntryPointNotFoundException)
                {
                }
                catch (BadImageFormatException)
                {
                }
            }
        }

        private static string GetPlatformString(IntPtr platform, uint param)
        {
            if (clGetPlatformInfo(platform, param, 0, null, out nuint size) != Success || size == 0)
                return "";
            var buffer = new byte[(int)size];
            clGetPlatformInfo(platform, param, size, buffer, out _);
            return Encoding.ASCII.GetString(buffer).TrimEnd('\0');
        }

        private static string GetDeviceString(IntPtr device, uint param)
        {
            if (clGetDeviceInfo(device, param, 0, null, out nuint size) != Success || size == 0)
                return "";
            var buffer = new byte[(int)size];
            clGetDeviceInfo(device, param, size, buffer, out _);
            return Encoding.ASCII.GetString(buffer).TrimEnd('\0');
        }

        private static ulong GetDeviceULong(IntPtr device, uint param)
        {
            var buffer = new byte[8];
            return clGetDeviceInfo(device, param, 8, buffer, out _) == Success ? BitConverter.ToUInt64(buffer, 0) : 0;
        }

        private static uint GetDeviceUInt(IntPtr device, uint param)
        {
            var buffer = new byte[4];
            return clGetDeviceInfo(device, param, 4, buffer, out _) == Success ? BitConverter.ToUInt32(buffer, 0) : 0;
        }

        private static ulong GetDeviceSize(IntPtr device, uint param)
        {
            var buffer = new byte[IntPtr.Size];
            if (clGetDeviceInfo(device, param, (nuint)IntPtr.Size, buffer, out _) != Success)
                return 0;
            return IntPtr.Size == 8 ? BitConverter.ToUInt64(buffer, 0) : BitConverter.ToUInt32(buffer, 0);
        }

        internal static void Check(int code, string call)
        {
            if (code != Success)
                throw new OpenCLException(call, code);
        }

        internal static string ErrorName(int code) => code switch
        {
            0 => "CL_SUCCESS",
            -1 => "CL_DEVICE_NOT_FOUND",
            -2 => "CL_DEVICE_NOT_AVAILABLE",
            -3 => "CL_COMPILER_NOT_AVAILABLE",
            -4 => "CL_MEM_OBJECT_ALLOCATION_FAILURE",
            -5 => "CL_OUT_OF_RESOURCES",
            -6 => "CL_OUT_OF_HOST_MEMORY",
            -11 => "CL_BUILD_PROGRAM_FAILURE",
            -30 => "CL_INVALID_VALUE",
            -33 => "CL_INVALID_DEVICE",
            -34 => "CL_INVALID_CONTEXT",
            -36 => "CL_INVALID_COMMAND_QUEUE",
            -38 => "CL_INVALID_MEM_OBJECT",
            -44 => "CL_INVALID_PROGRAM",
            -46 => "CL_INVALID_KERNEL_NAME",
            -48 => "CL_INVALID_KERNEL",
            -49 => "CL_INVALID_ARG_INDEX",
            -50 => "CL_INVALID_ARG_VALUE",
            -51 => "CL_INVALID_ARG_SIZE",
            -52 => "CL_INVALID_KERNEL_ARGS",
            -54 => "CL_INVALID_WORK_GROUP_SIZE",
            -61 => "CL_INVALID_BUFFER_SIZE",
            -63 => "CL_INVALID_GLOBAL_WORK_SIZE",
            -1001 => "CL_PLATFORM_NOT_FOUND_KHR",
            _ => "CL_ERROR"
        };

        #endregion

        #region Native Functions

        [DllImport(Library)] internal static extern int clGetPlatformIDs(uint numEntries, IntPtr[]? platforms, out uint numPlatforms);
        [DllImport(Library)] internal static extern int clGetPlatformInfo(IntPtr platform, uint param, nuint size, byte[]? value, out nuint sizeRet);
        [DllImport(Library)] internal static extern int clGetDeviceIDs(IntPtr platform, ulong type, uint numEntries, IntPtr[]? devices, out uint numDevices);
        [DllImport(Library)] internal static extern int clGetDeviceInfo(IntPtr device, uint param, nuint size, byte[]? value, out nuint sizeRet);

        [DllImport(Library)] internal static extern IntPtr clCreateContext(IntPtr[] properties, uint numDevices, IntPtr[] devices, IntPtr callback, IntPtr userData, out int error);
        [DllImport(Library)] internal static extern IntPtr clCreateCommandQueue(IntPtr context, IntPtr device, ulong properties, out int error);

        [DllImport(Library)]
        internal static extern IntPtr clCreateProgramWithSource(IntPtr context, uint count,
            [MarshalAs(UnmanagedType.LPArray, ArraySubType = UnmanagedType.LPStr)] string[] sources, IntPtr lengths, out int error);

        [DllImport(Library)]
        internal static extern int clBuildProgram(IntPtr program, uint numDevices, IntPtr[] devices,
            [MarshalAs(UnmanagedType.LPStr)] string options, IntPtr callback, IntPtr userData);

        [DllImport(Library)] internal static extern int clGetProgramBuildInfo(IntPtr program, IntPtr device, uint param, nuint size, byte[]? value, out nuint sizeRet);
        [DllImport(Library)] internal static extern IntPtr clCreateKernel(IntPtr program, [MarshalAs(UnmanagedType.LPStr)] string name, out int error);

        [DllImport(Library)] internal static extern int clSetKernelArg(IntPtr kernel, uint index, nuint size, ref IntPtr value);
        [DllImport(Library)] internal static extern int clSetKernelArg(IntPtr kernel, uint index, nuint size, ref int value);
        [DllImport(Library)] internal static extern int clSetKernelArg(IntPtr kernel, uint index, nuint size, ref float value);

        [DllImport(Library)] internal static extern IntPtr clCreateBuffer(IntPtr context, ulong flags, nuint size, IntPtr hostPtr, out int error);

        [DllImport(Library)] internal static extern int clEnqueueWriteBuffer(IntPtr queue, IntPtr buffer, uint blocking, nuint offset, nuint size, float[] data, uint numEvents, IntPtr waitList, IntPtr evt);
        [DllImport(Library)] internal static extern int clEnqueueWriteBuffer(IntPtr queue, IntPtr buffer, uint blocking, nuint offset, nuint size, int[] data, uint numEvents, IntPtr waitList, IntPtr evt);
        [DllImport(Library)] internal static extern int clEnqueueWriteBuffer(IntPtr queue, IntPtr buffer, uint blocking, nuint offset, nuint size, byte[] data, uint numEvents, IntPtr waitList, IntPtr evt);
        [DllImport(Library)] internal static extern int clEnqueueReadBuffer(IntPtr queue, IntPtr buffer, uint blocking, nuint offset, nuint size, float[] data, uint numEvents, IntPtr waitList, IntPtr evt);
        [DllImport(Library)] internal static extern int clEnqueueReadBuffer(IntPtr queue, IntPtr buffer, uint blocking, nuint offset, nuint size, int[] data, uint numEvents, IntPtr waitList, IntPtr evt);
        [DllImport(Library)] internal static extern int clEnqueueReadBuffer(IntPtr queue, IntPtr buffer, uint blocking, nuint offset, nuint size, byte[] data, uint numEvents, IntPtr waitList, IntPtr evt);

        [DllImport(Library)] internal static extern int clEnqueueFillBuffer(IntPtr queue, IntPtr buffer, ref int pattern, nuint patternSize, nuint offset, nuint size, uint numEvents, IntPtr waitList, IntPtr evt);
        [DllImport(Library)] internal static extern int clEnqueueCopyBuffer(IntPtr queue, IntPtr source, IntPtr destination, nuint sourceOffset, nuint destinationOffset, nuint size, uint numEvents, IntPtr waitList, IntPtr evt);

        [DllImport(Library)]
        internal static extern int clEnqueueNDRangeKernel(IntPtr queue, IntPtr kernel, uint workDim, IntPtr globalOffset,
            nuint[] globalSize, nuint[]? localSize, uint numEvents, IntPtr waitList, IntPtr evt);

        [DllImport(Library)] internal static extern int clFinish(IntPtr queue);
        [DllImport(Library)] internal static extern int clReleaseMemObject(IntPtr buffer);
        [DllImport(Library)] internal static extern int clReleaseKernel(IntPtr kernel);
        [DllImport(Library)] internal static extern int clReleaseProgram(IntPtr program);
        [DllImport(Library)] internal static extern int clReleaseCommandQueue(IntPtr queue);
        [DllImport(Library)] internal static extern int clReleaseContext(IntPtr context);

        #endregion
    }

    /// <summary>
    /// An OpenCL context, command queue and compiled program on one device, with
    /// small helpers for buffers and kernel launches.
    /// </summary>
    internal sealed class OpenCLRuntime : IDisposable
    {
        private readonly Dictionary<string, IntPtr> _kernels = new();
        private readonly List<IntPtr> _buffers = new();
        private IntPtr _context;
        private IntPtr _queue;
        private IntPtr _program;

        public OpenCLDeviceInfo Device { get; }

        public OpenCLRuntime(OpenCLDeviceInfo device, string source, string options = "")
        {
            Device = device;

            var properties = new[] { (IntPtr)0x1084 /* CL_CONTEXT_PLATFORM */, device.Platform, IntPtr.Zero };
            _context = OpenCL.clCreateContext(properties, 1, new[] { device.Device }, IntPtr.Zero, IntPtr.Zero, out int error);
            OpenCL.Check(error, "clCreateContext");

            _queue = OpenCL.clCreateCommandQueue(_context, device.Device, 0, out error);
            OpenCL.Check(error, "clCreateCommandQueue");

            _program = OpenCL.clCreateProgramWithSource(_context, 1, new[] { source }, IntPtr.Zero, out error);
            OpenCL.Check(error, "clCreateProgramWithSource");

            error = OpenCL.clBuildProgram(_program, 1, new[] { device.Device }, options, IntPtr.Zero, IntPtr.Zero);
            if (error != OpenCL.Success)
                throw new OpenCLException("clBuildProgram", error, BuildLog());
        }

        private string BuildLog()
        {
            if (OpenCL.clGetProgramBuildInfo(_program, Device.Device, OpenCL.ProgramBuildLog, 0, null, out nuint size) != OpenCL.Success || size == 0)
                return "";
            var buffer = new byte[(int)size];
            OpenCL.clGetProgramBuildInfo(_program, Device.Device, OpenCL.ProgramBuildLog, size, buffer, out _);
            return Encoding.UTF8.GetString(buffer).TrimEnd('\0');
        }

        public IntPtr Kernel(string name)
        {
            if (_kernels.TryGetValue(name, out var kernel))
                return kernel;

            kernel = OpenCL.clCreateKernel(_program, name, out int error);
            OpenCL.Check(error, "clCreateKernel(" + name + ")");
            _kernels[name] = kernel;
            return kernel;
        }

        public IntPtr CreateBuffer(long bytes)
        {
            var buffer = OpenCL.clCreateBuffer(_context, OpenCL.MemReadWrite, (nuint)Math.Max(bytes, 4), IntPtr.Zero, out int error);
            OpenCL.Check(error, "clCreateBuffer");
            _buffers.Add(buffer);
            return buffer;
        }

        public void ReleaseBuffer(ref IntPtr buffer)
        {
            if (buffer == IntPtr.Zero)
                return;
            _buffers.Remove(buffer);
            OpenCL.clReleaseMemObject(buffer);
            buffer = IntPtr.Zero;
        }

        public void Write(IntPtr buffer, float[] data, int count) => OpenCL.Check(
            count <= 0 ? 0 : OpenCL.clEnqueueWriteBuffer(_queue, buffer, 1, 0, (nuint)(count * 4L), data, 0, IntPtr.Zero, IntPtr.Zero), "clEnqueueWriteBuffer");

        public void Write(IntPtr buffer, int[] data, int count) => OpenCL.Check(
            count <= 0 ? 0 : OpenCL.clEnqueueWriteBuffer(_queue, buffer, 1, 0, (nuint)(count * 4L), data, 0, IntPtr.Zero, IntPtr.Zero), "clEnqueueWriteBuffer");

        public void Write(IntPtr buffer, byte[] data, int count) => OpenCL.Check(
            count <= 0 ? 0 : OpenCL.clEnqueueWriteBuffer(_queue, buffer, 1, 0, (nuint)count, data, 0, IntPtr.Zero, IntPtr.Zero), "clEnqueueWriteBuffer");

        public void Read(IntPtr buffer, float[] data, int count) => OpenCL.Check(
            count <= 0 ? 0 : OpenCL.clEnqueueReadBuffer(_queue, buffer, 1, 0, (nuint)(count * 4L), data, 0, IntPtr.Zero, IntPtr.Zero), "clEnqueueReadBuffer");

        public void Read(IntPtr buffer, int[] data, int count) => OpenCL.Check(
            count <= 0 ? 0 : OpenCL.clEnqueueReadBuffer(_queue, buffer, 1, 0, (nuint)(count * 4L), data, 0, IntPtr.Zero, IntPtr.Zero), "clEnqueueReadBuffer");

        public void Read(IntPtr buffer, byte[] data, int count) => OpenCL.Check(
            count <= 0 ? 0 : OpenCL.clEnqueueReadBuffer(_queue, buffer, 1, 0, (nuint)count, data, 0, IntPtr.Zero, IntPtr.Zero), "clEnqueueReadBuffer");

        public void FillInt(IntPtr buffer, int value, int count)
        {
            if (count <= 0)
                return;
            OpenCL.Check(OpenCL.clEnqueueFillBuffer(_queue, buffer, ref value, 4, 0, (nuint)(count * 4L), 0, IntPtr.Zero, IntPtr.Zero),
                "clEnqueueFillBuffer");
        }

        public void CopyInt(IntPtr source, IntPtr destination, int count)
        {
            if (count <= 0)
                return;
            OpenCL.Check(OpenCL.clEnqueueCopyBuffer(_queue, source, destination, 0, 0, (nuint)(count * 4L), 0, IntPtr.Zero, IntPtr.Zero),
                "clEnqueueCopyBuffer");
        }

        /// <summary>
        /// Sets the arguments of a kernel and launches it over <paramref name="globalSize"/>
        /// work items (rounded up to the work-group size). Arguments may be buffers
        /// (<see cref="IntPtr"/>), <see cref="int"/> or <see cref="float"/> values.
        /// </summary>
        public void Run(string kernelName, int globalSize, int localSize, params object[] arguments)
        {
            if (globalSize <= 0)
                return;

            var kernel = Kernel(kernelName);
            for (uint i = 0; i < arguments.Length; i++)
            {
                int error;
                switch (arguments[i])
                {
                    case IntPtr buffer:
                        error = OpenCL.clSetKernelArg(kernel, i, (nuint)IntPtr.Size, ref buffer);
                        break;
                    case int value:
                        error = OpenCL.clSetKernelArg(kernel, i, 4, ref value);
                        break;
                    case float value:
                        error = OpenCL.clSetKernelArg(kernel, i, 4, ref value);
                        break;
                    default:
                        throw new ArgumentException($"Unsupported kernel argument type {arguments[i]?.GetType().Name}.");
                }
                OpenCL.Check(error, $"clSetKernelArg({kernelName}, {i})");
            }

            long rounded = (globalSize + (long)localSize - 1) / localSize * localSize;
            OpenCL.Check(OpenCL.clEnqueueNDRangeKernel(_queue, kernel, 1, IntPtr.Zero,
                new[] { (nuint)rounded }, new[] { (nuint)localSize }, 0, IntPtr.Zero, IntPtr.Zero),
                "clEnqueueNDRangeKernel(" + kernelName + ")");
        }

        public void Finish() => OpenCL.Check(OpenCL.clFinish(_queue), "clFinish");

        public void Dispose()
        {
            foreach (var buffer in _buffers)
                OpenCL.clReleaseMemObject(buffer);
            _buffers.Clear();

            foreach (var kernel in _kernels.Values)
                OpenCL.clReleaseKernel(kernel);
            _kernels.Clear();

            if (_program != IntPtr.Zero) { OpenCL.clReleaseProgram(_program); _program = IntPtr.Zero; }
            if (_queue != IntPtr.Zero) { OpenCL.clReleaseCommandQueue(_queue); _queue = IntPtr.Zero; }
            if (_context != IntPtr.Zero) { OpenCL.clReleaseContext(_context); _context = IntPtr.Zero; }
        }
    }
}
#endif
