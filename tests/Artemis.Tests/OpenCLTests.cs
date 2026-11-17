using System.Numerics;
using Artemis.Compute;
using Artemis.Core;
using Artemis.Particles;

namespace Artemis.Tests;

/// <summary>
/// A test that needs an OpenCL device (a GPU, or a CPU runtime such as POCL). It is
/// skipped on machines without one.
/// </summary>
public sealed class OpenCLFactAttribute : FactAttribute
{
    public OpenCLFactAttribute()
    {
        if (OpenCL.SelectDevice(allowCpu: true) == null)
            Skip = "No OpenCL device available (install a GPU driver, or pocl-opencl-icd for a CPU device).";
    }
}

public class OpenCLTests
{
    private const float Dt = 1f / 60f;

    private static MassiveParticleSystem Create(bool gpu, int capacity = 20_000)
    {
        var system = new MassiveParticleSystem(capacity, new AABB(new Vector3D(-10, 0, -10), new Vector3D(10, 20, 10)))
        {
            Friction = 0.5f,
            Restitution = 0.1f,
            AllowOpenCLCpuDevice = true
        };
        return system;
    }

    private static void EnableGpu(MassiveParticleSystem system)
    {
        system.UseGPU = true;
        Assert.True(system.UseGPU, system.GpuError);
    }

    private static void AddColumn(MassiveParticleSystem system, int nx, int ny, int nz, float r = 0.05f)
    {
        var rng = new Random(3);
        for (int y = 0; y < ny; y++)
            for (int z = 0; z < nz; z++)
                for (int x = 0; x < nx; x++)
                    system.AddParticle(new Vector3(
                        (x - nx / 2) * 2 * r + ((float)rng.NextDouble() - 0.5f) * 0.03f,
                        r + y * 2 * r + ((float)rng.NextDouble() - 0.5f) * 0.03f,
                        (z - nz / 2) * 2 * r + ((float)rng.NextDouble() - 0.5f) * 0.03f),
                        Vector3.Zero, 1f, r);
    }

    private static void Run(MassiveParticleSystem system, double seconds)
    {
        for (int i = 0; i < (int)Math.Round(seconds / Dt); i++)
            system.Update(Dt);
    }

    [Fact]
    public void Without_a_device_the_gpu_switch_stays_off_and_says_why()
    {
        var system = new MassiveParticleSystem(100) { AllowOpenCLCpuDevice = false };
        system.UseGPU = true;

        if (OpenCL.SelectDevice(allowCpu: false) == null)
        {
            Assert.False(system.UseGPU);
            Assert.False(string.IsNullOrEmpty(system.GpuError));
            Assert.StartsWith("CPU", system.GpuBackendInfo);
        }
        else
        {
            Assert.True(system.UseGPU);
            Assert.StartsWith("OpenCL", system.GpuBackendInfo);
        }
    }

    [OpenCLFact]
    public void Kernels_compile_and_the_first_step_matches_the_cpu_solver()
    {
        var cpu = Create(false);
        var gpu = Create(true);
        AddColumn(cpu, 12, 40, 12);
        AddColumn(gpu, 12, 40, 12);
        EnableGpu(gpu);

        cpu.Update(Dt);
        gpu.Update(Dt);

        Assert.True(gpu.LastUpdateUsedGpu);
        float worst = 0;
        for (int i = 0; i < cpu.Count; i++)
            worst = MathF.Max(worst, Vector3.Distance(cpu.GetParticle(i).Position, gpu.GetParticle(i).Position));
        Assert.True(worst < 1e-4f, $"GPU and CPU differ by {worst} m after one frame");
        Assert.Equal(cpu.ActiveCount, gpu.ActiveCount);
        Assert.Equal(cpu.AliveCount, gpu.AliveCount);
    }

    [OpenCLFact]
    public void Sand_column_collapses_into_the_same_kind_of_heap_on_the_gpu()
    {
        var gpu = Create(true);
        AddColumn(gpu, 10, 60, 10);
        EnableGpu(gpu);

        Run(gpu, 25);

        float top = 0;
        var radii = new List<float>();
        for (int i = 0; i < gpu.Count; i++)
        {
            var p = gpu.GetParticle(i).Position;
            Assert.False(float.IsNaN(p.X));
            Assert.True(p.Y >= 0.049f);
            top = MathF.Max(top, p.Y);
            radii.Add(MathF.Sqrt(p.X * p.X + p.Z * p.Z));
        }
        radii.Sort();
        float reach = radii[(int)(radii.Count * 0.98)];

        Assert.InRange(top, 0.5f, 3.0f);
        Assert.InRange(Math.Atan2(top, reach) * 180 / Math.PI, 15, 45);
        Assert.True(gpu.SleepingCount > gpu.Count * 0.95, $"only {gpu.SleepingCount} grains asleep");
    }

    [OpenCLFact]
    public void Gpu_simulation_is_deterministic()
    {
        Vector3[] Simulate()
        {
            var system = Create(true);
            AddColumn(system, 10, 30, 10);
            EnableGpu(system);
            Run(system, 2);
            return Enumerable.Range(0, system.Count).Select(i => system.GetParticle(i).Position).ToArray();
        }

        var a = Simulate();
        var b = Simulate();
        for (int i = 0; i < a.Length; i++)
            Assert.Equal(a[i], b[i]);
    }

    [OpenCLFact]
    public void Host_edits_reach_the_device_and_device_results_reach_the_host()
    {
        var system = Create(true);
        var statue = TriangleMesh.CreateBox(new Vector3(-0.5f, 0, -0.5f), new Vector3(0.5f, 3, 0.5f));
        int count = system.SpawnFromMesh(statue, 0.05f, out _);
        EnableGpu(system);

        Run(system, 0.5);
        Assert.Equal(count, system.FrozenCount);        // read back from the device
        Assert.Equal(0, system.ActiveCount);

        system.ReleaseBeyondPlane(Vector3.UnitY, 2.0f);   // host edit, uploaded on the next update
        Run(system, 0.5);
        Assert.True(system.ActiveCount > 0);
        Assert.True(system.FrozenCount < count);

        system.ReleaseAll();
        Run(system, 8);
        float top = 0;
        for (int i = 0; i < count; i++) top = MathF.Max(top, system.GetParticle(i).Position.Y);
        Assert.True(top < 2.0f, $"released statue still {top} m tall");
        Assert.True(system.LastUpdateUsedGpu);
    }

    [OpenCLFact]
    public void Lifetimes_expire_on_the_gpu()
    {
        var system = Create(true);
        for (int i = 0; i < 100; i++)
            system.AddParticle(new Vector3(i * 0.2f - 10, 5, 0), Vector3.Zero, 1f, 0.05f, lifetime: i < 50 ? 0.25f : float.MaxValue);
        EnableGpu(system);

        Run(system, 0.5);

        Assert.Equal(50, system.AliveCount);
        Assert.False(system.GetParticle(0).IsAlive);
        Assert.True(system.GetParticle(99).IsAlive);
    }

    [OpenCLFact]
    public void Soft_bodies_run_on_the_gpu_and_match_the_cpu()
    {
        Vector3 Simulate(bool gpu, out bool usedGpu)
        {
            var system = Create(gpu);
            system.SpawnSoftBodyFromMesh(TriangleMesh.CreateBox(new Vector3(-0.5f, 1, -0.5f), new Vector3(0.5f, 2, 0.5f)),
                0.05f, out int body, stiffness: 0.3f);
            if (gpu) EnableGpu(system);
            Run(system, 3);
            usedGpu = system.LastUpdateUsedGpu;
            return system.GetSoftBodyCenter(body);
        }

        var onCpu = Simulate(false, out _);
        var onGpu = Simulate(true, out bool used);

        Assert.True(used);
        Assert.True(Vector3.Distance(onCpu, onGpu) < 0.005f, $"cpu {onCpu} vs gpu {onGpu}");
    }

    [OpenCLFact]
    public void Rigid_groups_fall_back_to_the_cpu_and_still_shatter()
    {
        var system = Create(true, 50_000);
        var profile = new List<Vector2> { new(0.025f, 0), new(0.03f, 0.02f), new(0.035f, 0.10f) };
        var glass = TriangleMesh.CreateLathe(profile, 48, 0.004f).Transform(1f, new Vector3(0, 1, 0));
        system.SpawnRigidBodyFromMesh(glass, 0.002f, 2.0f, out int id, mass: 0.001f);
        EnableGpu(system);

        Run(system, 1.5);

        Assert.False(system.LastUpdateUsedGpu);
        Assert.False(system.IsGroupAlive(id));
        Assert.True(system.RigidGroupCount >= 10);
    }
}
