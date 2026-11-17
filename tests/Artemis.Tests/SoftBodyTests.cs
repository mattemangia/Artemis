using System.Numerics;
using Artemis.Core;
using Artemis.Particles;

namespace Artemis.Tests;

public class SoftBodyTests
{
    private const float Dt = 1f / 60f;

    private static MassiveParticleSystem World(float gravity = -9.81f)
        => new(40_000, new AABB(new Vector3D(-5, 0, -5), new Vector3D(5, 10, 5)))
        {
            Friction = 0.5f,
            Restitution = 0.1f,
            Gravity = new Vector3(0, gravity, 0)
        };

    private static TriangleMesh Cube(float size, float bottom)
        => TriangleMesh.CreateBox(new Vector3(-size / 2, bottom, -size / 2), new Vector3(size / 2, bottom + size, size / 2));

    private static (float Min, float Max) VerticalExtent(MassiveParticleSystem system, int first, int count)
    {
        float min = float.MaxValue, max = float.MinValue;
        for (int i = first; i < first + count; i++)
        {
            float y = system.GetParticle(i).Position.Y;
            min = MathF.Min(min, y);
            max = MathF.Max(max, y);
        }
        return (min, max);
    }

    private static void Run(MassiveParticleSystem system, double seconds)
    {
        for (int i = 0; i < (int)Math.Round(seconds / Dt); i++)
            system.Update(Dt);
    }

    [Fact]
    public void Soft_cube_squashes_on_impact_and_springs_back_to_its_shape()
    {
        var system = World();
        int n = system.SpawnSoftBodyFromMesh(Cube(1f, 1f), 0.05f, out int body, stiffness: 0.3f);
        Assert.Equal(1000, n);
        Assert.True(system.GetSoftBodyClusterCount(body) > 50);

        float squashed = float.MaxValue;
        for (int f = 0; f < 360; f++)
        {
            system.Update(Dt);
            var (min, max) = VerticalExtent(system, 0, n);
            if (f > 20) squashed = MathF.Min(squashed, max - min);
        }

        var (finalMin, finalMax) = VerticalExtent(system, 0, n);
        Assert.True(squashed < 0.8f, $"the cube barely deformed (min height {squashed})");
        Assert.InRange(finalMax - finalMin, 0.82f, 0.92f);    // 0.9 tall at rest, sags a little
        Assert.InRange(finalMin, 0.049f, 0.06f);               // resting on the floor
        Assert.True(system.GetSoftBodyVelocity(body).Length() < 0.02f, "the cube is still moving");

        // Symmetric body, symmetric drop: it must not drift sideways.
        var center = system.GetSoftBodyCenter(body);
        Assert.InRange(center.X, -0.01f, 0.01f);
        Assert.InRange(center.Z, -0.01f, 0.01f);
    }

    [Fact]
    public void Stiffer_bodies_sag_less_under_their_own_weight()
    {
        float Height(float stiffness)
        {
            var system = World();
            int n = system.SpawnSoftBodyFromMesh(Cube(1f, 0f), 0.05f, out _, stiffness: stiffness, damping: 0.5f);
            Run(system, 4);
            var (min, max) = VerticalExtent(system, 0, n);
            return max - min;
        }

        float jelly = Height(0.05f), rubber = Height(0.3f), hard = Height(0.95f);
        Assert.True(jelly < rubber && rubber < hard, $"heights {jelly}, {rubber}, {hard}");
        Assert.InRange(hard, 0.88f, 0.905f);
        Assert.True(jelly > 0.5f, "jelly collapsed");
    }

    [Fact]
    public void Damping_reduces_the_rebound_without_slowing_the_fall()
    {
        (float Bounce, float LandingTime) Drop(float damping)
        {
            var system = World();
            system.SpawnSoftBodyFromMesh(Cube(1f, 1f), 0.05f, out int body, stiffness: 0.3f, damping: damping);
            float peak = 0, landing = -1;
            bool rising = false;
            for (int f = 0; f < 120; f++)
            {
                system.Update(Dt);
                var v = system.GetSoftBodyVelocity(body);
                if (landing < 0 && v.Y > -0.5f && f > 10) landing = f * Dt;
                if (v.Y > 0) rising = true;
                if (rising) peak = MathF.Max(peak, system.GetSoftBodyCenter(body).Y);
            }
            return (peak, landing);
        }

        var rubber = Drop(0f);
        var dead = Drop(0.5f);

        Assert.True(rubber.Bounce > 0.75f, $"an undamped rubber cube should bounce (peak {rubber.Bounce})");
        Assert.True(dead.Bounce < 0.55f, $"a damped cube should barely bounce (peak {dead.Bounce})");
        // Free fall is not affected by damping: both land at the same moment.
        Assert.InRange(dead.LandingTime - rubber.LandingTime, -0.05f, 0.05f);
    }

    [Fact]
    public void Damping_does_not_slow_down_linear_or_angular_motion()
    {
        var system = new MassiveParticleSystem(10_000, new AABB(new Vector3D(-50, -50, -50), new Vector3D(50, 50, 50)))
        {
            Gravity = Vector3.Zero,
            Damping = 1f
        };
        var bar = TriangleMesh.CreateBox(new Vector3(-0.6f, -0.1f, -0.1f), new Vector3(0.6f, 0.1f, 0.1f));
        system.SpawnSoftBodyFromMesh(bar, 0.025f, out int body, stiffness: 0.5f, damping: 0.5f);
        system.SetSoftBodyVelocity(body, new Vector3(0.5f, 0, 0), new Vector3(0, 0, 2f));

        Run(system, 2);

        var center = system.GetSoftBodyCenter(body);
        Assert.InRange(center.X, 0.97f, 1.03f);
        Assert.InRange(system.GetSoftBodyVelocity(body).X, 0.49f, 0.51f);
        // The spinning bar stretches a little, which slows it a little; it must not stop.
        Assert.InRange(system.GetSoftBodyAngularVelocity(body).Z, 1.8f, 2.05f);
    }

    [Fact]
    public void Soft_bodies_stack_without_sinking_into_each_other()
    {
        var system = World();
        int a = system.SpawnSoftBodyFromMesh(Cube(1f, 0f), 0.05f, out int lower, stiffness: 0.5f);
        int b = system.SpawnSoftBodyFromMesh(Cube(0.6f, 1.2f), 0.05f, out int upper, stiffness: 0.5f);

        Run(system, 5);

        var (lowerMin, lowerMax) = VerticalExtent(system, 0, a);
        var (upperMin, _) = VerticalExtent(system, a, b);

        // Particle centres of touching bodies are one diameter (0.1) apart.
        Assert.InRange(upperMin - lowerMax, 0.06f, 0.12f);
        Assert.InRange(lowerMin, 0.049f, 0.06f);
        Assert.True(system.GetSoftBodyVelocity(upper).Length() < 0.02f);
    }

    [Fact]
    public void Stiffness_and_damping_can_be_changed_while_running()
    {
        var system = World();
        int n = system.SpawnSoftBodyFromMesh(Cube(1f, 0f), 0.05f, out int body, stiffness: 0.9f, damping: 0.5f);
        Run(system, 2);
        var (min0, max0) = VerticalExtent(system, 0, n);

        system.SetSoftBodyStiffness(body, 0.02f);
        Assert.Equal(0.02f, system.GetSoftBodyStiffness(body));
        Run(system, 3);
        var (min1, max1) = VerticalExtent(system, 0, n);

        Assert.True(max1 - min1 < max0 - min0 - 0.05f, "going limp should make the body sag");
    }

    [Fact]
    public void Dissolving_a_soft_body_turns_it_into_loose_grains()
    {
        var system = World();
        int n = system.SpawnSoftBodyFromMesh(Cube(1f, 0f), 0.05f, out int body, stiffness: 0.5f);
        Run(system, 1);

        system.DissolveSoftBody(body);
        Assert.Equal(0, system.SoftBodyCount);
        Assert.False(system.IsSoftBodyAlive(body));
        Assert.Equal(-1, system.GetParticleSoftBody(0));

        // Loose grains fall asleep once they rest; members of a soft body never do.
        Run(system, 4);
        Assert.True(system.SleepingCount > n * 0.9, $"only {system.SleepingCount} of {n} grains asleep");

        // And nothing holds them together any more.
        system.ApplyExplosion(new Vector3(0, 0.5f, 0), force: 8f, radius: 2f);
        Run(system, 2);
        var (_, max) = VerticalExtent(system, 0, n);
        Assert.True(max < 0.7f, $"the grains should have been scattered (top at {max})");
    }

    [Fact]
    public void Compact_keeps_soft_bodies_working()
    {
        var system = World();
        system.AddParticle(new Vector3(3, 1, 3), Vector3.Zero, 1f, 0.05f, lifetime: 0.1f);
        int n = system.SpawnSoftBodyFromMesh(Cube(1f, 1f), 0.05f, out int body, stiffness: 0.5f);
        Run(system, 0.5);

        system.Compact(); // removes the expired particle in front of the body
        Assert.Equal(n, system.Count);
        Assert.Equal(body, system.GetParticleSoftBody(0));

        Run(system, 4);
        var (min, max) = VerticalExtent(system, 0, n);
        Assert.InRange(max - min, 0.85f, 0.92f);
        Assert.InRange(min, 0.049f, 0.06f);
    }

    [Fact]
    public void Soft_body_simulation_is_deterministic()
    {
        static Vector3[] Simulate()
        {
            var system = World();
            system.SpawnSoftBodyFromMesh(Cube(1f, 1f), 0.05f, out _, stiffness: 0.2f);
            system.SpawnSoftBodyFromMesh(Cube(0.5f, 2.5f), 0.05f, out _, stiffness: 0.6f);
            Run(system, 2);
            var result = new Vector3[system.Count];
            for (int i = 0; i < result.Length; i++) result[i] = system.GetParticle(i).Position;
            return result;
        }

        var a = Simulate();
        var b = Simulate();
        for (int i = 0; i < a.Length; i++)
            Assert.Equal(a[i], b[i]);
    }
}
