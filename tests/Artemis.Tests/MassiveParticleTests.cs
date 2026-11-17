using System.Numerics;
using Artemis.Core;
using Artemis.Particles;

namespace Artemis.Tests;

public class MassiveParticleTests
{
    private const float Dt = 1f / 60f;

    private static MassiveParticleSystem Create(int capacity, float halfWidth = 10, float height = 20)
        => new(capacity, new AABB(new Vector3D(-halfWidth, 0, -halfWidth), new Vector3D(halfWidth, height, halfWidth)))
        {
            Friction = 0.5f,
            Restitution = 0.1f
        };

    private static void Run(MassiveParticleSystem system, double seconds)
    {
        for (int i = 0; i < (int)Math.Round(seconds / Dt); i++)
            system.Update(Dt);
    }

    private static void AddBlock(MassiveParticleSystem system, int nx, int ny, int nz, float r, float jitter, float baseY = 0)
    {
        var rng = new Random(3);
        for (int y = 0; y < ny; y++)
            for (int z = 0; z < nz; z++)
                for (int x = 0; x < nx; x++)
                {
                    system.AddParticle(new Vector3(
                        (x - nx / 2) * 2 * r + ((float)rng.NextDouble() - 0.5f) * jitter,
                        baseY + r + y * 2 * r + ((float)rng.NextDouble() - 0.5f) * jitter,
                        (z - nz / 2) * 2 * r + ((float)rng.NextDouble() - 0.5f) * jitter),
                        Vector3.Zero, 1f, r);
                }
    }

    [Fact]
    public void Free_falling_particle_follows_the_analytic_trajectory()
    {
        var system = new MassiveParticleSystem(16, new AABB(new Vector3D(-10, -1000, -10), new Vector3D(10, 1000, 10)))
        {
            Damping = 1f
        };
        system.AddParticle(new Vector3(0, 100, 0), new Vector3(2, 0, 0), 1f, 0.05f);

        Run(system, 2);

        var p = system.GetParticle(0);
        Assert.InRange(100 - p.Position.Y, 19.5f, 19.8f); // 0.5 g t^2 = 19.62
        Assert.Equal(4.0f, p.Position.X, 2);
        Assert.Equal(-19.62f, p.Velocity.Y, 1);
    }

    [Fact]
    public void Intersecting_spawn_is_separated_without_blowing_apart()
    {
        var system = Create(4096);
        // 2000 grains thrown into a volume a third of what they need.
        var rng = new Random(5);
        for (int i = 0; i < 2000; i++)
        {
            system.AddParticle(new Vector3((float)rng.NextDouble() * 0.8f - 0.4f, 0.05f + (float)rng.NextDouble() * 0.8f,
                (float)rng.NextDouble() * 0.8f - 0.4f), Vector3.Zero, 1f, 0.05f);
        }

        float maxSpeed = 0;
        for (int f = 0; f < 240; f++)
        {
            system.Update(Dt);
            for (int i = 0; i < system.Count; i++)
                maxSpeed = MathF.Max(maxSpeed, system.GetParticle(i).Velocity.Length());
        }

        // Nothing may be launched: the fastest a grain can legitimately go is a
        // fall from the top of the heap.
        Assert.True(maxSpeed < 5f, $"a grain reached {maxSpeed} m/s");
        for (int i = 0; i < system.Count; i++)
        {
            var p = system.GetParticle(i).Position;
            Assert.False(float.IsNaN(p.X));
            Assert.InRange(p.Y, 0.049f, 3f);
        }
    }

    [Fact]
    public void Column_of_sand_collapses_into_a_heap_with_a_sand_like_slope_and_falls_asleep()
    {
        var system = Create(8192);
        AddBlock(system, 10, 60, 10, 0.05f, 0.03f); // 1 m wide, 6 m tall

        Run(system, 25);

        float top = 0, reach = 0;
        // Ignore the few grains that bounced away: measure where 98% of the sand is.
        var radii = new List<float>();
        for (int i = 0; i < system.Count; i++)
        {
            var p = system.GetParticle(i).Position;
            top = MathF.Max(top, p.Y);
            radii.Add(MathF.Sqrt(p.X * p.X + p.Z * p.Z));
        }
        radii.Sort();
        reach = radii[(int)(radii.Count * 0.98)];

        Assert.True(top < 3.0f, $"the column is still standing ({top} m tall)");
        Assert.True(top > 0.5f, $"the sand spread into a puddle ({top} m tall): no friction?");

        double slope = Math.Atan2(top, reach) * 180 / Math.PI;
        Assert.InRange(slope, 15, 45);

        Assert.True(system.SleepingCount > system.Count * 0.95, $"only {system.SleepingCount} of {system.Count} grains are asleep");
    }

    [Fact]
    public void Settled_heap_keeps_its_volume()
    {
        var system = Create(8192);
        AddBlock(system, 12, 12, 12, 0.05f, 0.02f);
        Run(system, 12);

        // Grains must not be squeezed into each other by the weight above them.
        int deepOverlaps = 0;
        var positions = new Vector3[system.Count];
        for (int i = 0; i < positions.Length; i++) positions[i] = system.GetParticle(i).Position;
        for (int i = 0; i < positions.Length; i++)
        {
            for (int j = i + 1; j < positions.Length; j++)
            {
                if (Vector3.DistanceSquared(positions[i], positions[j]) < 0.08f * 0.08f) // 20% overlap
                    deepOverlaps++;
            }
        }

        Assert.True(deepOverlaps < positions.Length / 100, $"{deepOverlaps} pairs overlap by more than 20%");
    }

    [Fact]
    public void Simulation_is_deterministic_across_runs()
    {
        static Vector3[] Simulate()
        {
            var system = Create(8192);
            AddBlock(system, 14, 30, 14, 0.05f, 0.03f);
            Run(system, 3);
            var result = new Vector3[system.Count];
            for (int i = 0; i < result.Length; i++) result[i] = system.GetParticle(i).Position;
            return result;
        }

        var a = Simulate();
        var b = Simulate();
        for (int i = 0; i < a.Length; i++)
            Assert.Equal(a[i], b[i]);
    }

    [Fact]
    public void Sleeping_sand_is_woken_when_its_support_is_removed()
    {
        var system = Create(8192);
        AddBlock(system, 6, 6, 6, 0.05f, 0.02f);
        system.Freeze(0, 36); // bottom layer acts as a shelf
        Run(system, 8);
        Assert.True(system.SleepingCount > 150, "the pile should be asleep on its shelf");

        // Lift the shelf away by moving the pile's world: release the shelf and
        // open the floor beneath it.
        system.Bounds = new AABB(new Vector3D(-10, -5, -10), new Vector3D(10, 20, 10));
        system.ReleaseAll();
        Run(system, 4);

        float highest = float.MinValue;
        for (int i = 0; i < system.Count; i++)
            highest = MathF.Max(highest, system.GetParticle(i).Position.Y);
        Assert.True(highest < -3f, $"grains were left hanging at y={highest}");
    }
}

public class MeshToParticlesTests
{
    private const float Dt = 1f / 60f;

    [Fact]
    public void Filling_a_box_gives_one_point_per_cell_of_its_volume()
    {
        var mesh = TriangleMesh.CreateBox(new Vector3(-1, 0, -0.5f), new Vector3(1, 3, 0.5f));
        var points = MeshSampler.FillVolume(mesh, 0.1f);

        Assert.Equal(20 * 30 * 10, points.Count);
        Assert.All(points, p =>
        {
            Assert.InRange(p.X, -1f, 1f);
            Assert.InRange(p.Y, 0f, 3f);
            Assert.InRange(p.Z, -0.5f, 0.5f);
        });
    }

    [Fact]
    public void Hollow_lathe_is_filled_only_inside_its_walls()
    {
        // A cylinder-shaped cup: radius 0.5, height 1, walls 0.1 thick.
        var profile = new List<Vector2> { new(0.5f, 0), new(0.5f, 1f) };
        var mesh = TriangleMesh.CreateLathe(profile, 64, wallThickness: 0.1f);
        var points = MeshSampler.FillVolume(mesh, 0.02f);

        Assert.NotEmpty(points);
        foreach (var p in points)
        {
            float radius = MathF.Sqrt(p.X * p.X + p.Z * p.Z);
            Assert.True(radius <= 0.5f, "point outside the cup");
            // Above the base, nothing may sit in the cavity.
            if (p.Y > 0.12f)
                Assert.True(radius >= 0.385f, $"point inside the cavity at radius {radius}, height {p.Y}");
        }

        // Volume: base disc + wall ring, to within the resolution of the 64-gon.
        double expected = Math.PI * 0.25 * 0.1 + Math.PI * (0.25 - 0.16) * 0.9;
        double actual = points.Count * 0.02 * 0.02 * 0.02;
        Assert.InRange(actual / expected, 0.93, 1.03);
    }

    [Fact]
    public void Obj_text_is_parsed_into_triangles()
    {
        const string obj = """
            # a quad and a triangle
            v 0 0 0
            v 1 0 0
            v 1 1 0
            v 0 1 0
            v 0.5 2 0
            f 1/1/1 2/2/1 3/3/1 4/4/1
            f -2 -3 -1
            """;

        var mesh = TriangleMesh.ParseObj(obj);

        Assert.Equal(5, mesh.Vertices.Count);
        Assert.Equal(3, mesh.TriangleCount); // quad -> 2 triangles, plus 1
        Assert.Equal(new[] { 0, 1, 2, 0, 2, 3, 3, 2, 4 }, mesh.Indices);
    }

    [Fact]
    public void Surface_sampling_covers_an_open_mesh()
    {
        var mesh = new TriangleMesh();
        mesh.Vertices.AddRange(new[] { new Vector3(0, 0, 0), new Vector3(1, 0, 0), new Vector3(1, 0, 1), new Vector3(0, 0, 1) });
        mesh.Indices.AddRange(new[] { 0, 1, 2, 0, 2, 3 });

        var points = MeshSampler.SampleSurface(mesh, 0.1f);

        Assert.InRange(points.Count, 100, 125); // a 1 m^2 sheet at 0.1 m spacing
    }

    [Fact]
    public void Mesh_made_of_frozen_grains_holds_its_shape_until_released_then_becomes_a_heap()
    {
        var system = new MassiveParticleSystem(20_000, new AABB(new Vector3D(-10, 0, -10), new Vector3D(10, 20, 10)))
        {
            Friction = 0.5f
        };
        var statue = TriangleMesh.CreateBox(new Vector3(-0.5f, 0, -0.5f), new Vector3(0.5f, 5, 0.5f));

        int count = system.SpawnFromMesh(statue, 0.05f, out int first);
        Assert.Equal(10 * 50 * 10, count);
        Assert.Equal(0, first);
        Assert.Equal(count, system.FrozenCount);

        for (int i = 0; i < 60; i++) system.Update(Dt);
        Assert.Equal(0, system.ActiveCount);
        float top = 0;
        for (int i = 0; i < count; i++) top = MathF.Max(top, system.GetParticle(i).Position.Y);
        Assert.InRange(top, 4.9f, 5.0f);

        Assert.Equal(count, system.ReleaseAll());
        for (int i = 0; i < 60 * 20; i++) system.Update(Dt);

        top = 0;
        for (int i = 0; i < count; i++) top = MathF.Max(top, system.GetParticle(i).Position.Y);
        Assert.True(top < 2.5f, $"the released statue is still {top} m tall");
        Assert.True(top > 0.4f, "the sand should pile up, not spread flat");
        Assert.Equal(0, system.FrozenCount);
    }

    [Fact]
    public void Statue_can_be_released_progressively_from_the_top()
    {
        var system = new MassiveParticleSystem(20_000, new AABB(new Vector3D(-10, 0, -10), new Vector3D(10, 20, 10)));
        var statue = TriangleMesh.CreateBox(new Vector3(-0.5f, 0, -0.5f), new Vector3(0.5f, 4, 0.5f));
        int count = system.SpawnFromMesh(statue, 0.05f, out _);

        int released = system.ReleaseBeyondPlane(Vector3.UnitY, 3.0f);

        Assert.InRange(released, count / 4 - 150, count / 4 + 150);
        Assert.Equal(count - released, system.FrozenCount);
    }
}

public class ShatteringBodyTests
{
    private const float Dt = 1f / 60f;

    private static MassiveParticleSystem Create()
        => new(50_000, new AABB(new Vector3D(-2, 0, -2), new Vector3D(2, 5, 2))) { Friction = 0.4f, Restitution = 0.1f };

    private static TriangleMesh Glass(float height)
    {
        var profile = new List<Vector2> { new(0.025f, 0), new(0.03f, 0.02f), new(0.035f, 0.10f) };
        return TriangleMesh.CreateLathe(profile, 48, 0.004f).Transform(1f, new Vector3(0, height, 0));
    }

    private static float MaxPairDistanceError(MassiveParticleSystem system, int[] indices, Vector3[] original)
    {
        float worst = 0;
        for (int a = 0; a < indices.Length; a += 7)
        {
            for (int b = a + 1; b < indices.Length; b += 13)
            {
                float now = Vector3.Distance(system.GetParticle(indices[a]).Position, system.GetParticle(indices[b]).Position);
                float before = Vector3.Distance(original[a], original[b]);
                worst = MathF.Max(worst, MathF.Abs(now - before));
            }
        }
        return worst;
    }

    [Fact]
    public void Rigid_body_falls_as_one_piece_and_rests_on_the_floor_without_sinking()
    {
        var system = Create();
        var box = TriangleMesh.CreateBox(new Vector3(-0.1f, 1, -0.1f), new Vector3(0.1f, 1.2f, 0.1f));
        int count = system.SpawnRigidBodyFromMesh(box, 0.01f, float.PositiveInfinity, out int group);
        Assert.Equal(1000, count);

        var indices = Enumerable.Range(0, count).ToArray();
        var original = indices.Select(i => system.GetParticle(i).Position).ToArray();

        for (int i = 0; i < 180; i++) system.Update(Dt);

        Assert.True(system.IsGroupAlive(group));
        Assert.True(MaxPairDistanceError(system, indices, original) < 1e-3f, "the body deformed");

        float lowest = float.MaxValue, highest = 0;
        for (int i = 0; i < count; i++)
        {
            float y = system.GetParticle(i).Position.Y;
            lowest = MathF.Min(lowest, y);
            highest = MathF.Max(highest, y);
        }
        Assert.InRange(lowest, 0.0099f, 0.0105f);   // bottom particles sit on the floor
        Assert.InRange(highest, 0.185f, 0.195f);    // and the box is still 0.2 m tall
    }

    [Fact]
    public void Spinning_rigid_body_keeps_spinning_without_gravity()
    {
        var system = new MassiveParticleSystem(5000, new AABB(new Vector3D(-50, -50, -50), new Vector3D(50, 50, 50)))
        {
            Gravity = Vector3.Zero,
            Damping = 1f
        };
        var bar = TriangleMesh.CreateBox(new Vector3(-0.3f, -0.05f, -0.05f), new Vector3(0.3f, 0.05f, 0.05f));
        int count = system.SpawnRigidBodyFromMesh(bar, 0.01f, float.PositiveInfinity, out int group);
        system.SetGroupVelocity(group, new Vector3(1, 0, 0), new Vector3(0, 0, 3));

        for (int i = 0; i < 120; i++) system.Update(Dt); // 2 s

        var (center, rotation) = system.GetGroupPose(group);
        Assert.InRange(center.X, 1.95f, 2.05f);

        // 3 rad/s for 2 s about Z.
        var axisX = Vector3.Transform(Vector3.UnitX, rotation);
        double angle = Math.Atan2(axisX.Y, axisX.X);
        double expected = Math.IEEERemainder(6.0, 2 * Math.PI);
        Assert.InRange(Math.Abs(Math.IEEERemainder(angle - expected, 2 * Math.PI)), 0, 0.15);
    }

    [Fact]
    public void Glass_dropped_from_one_metre_shatters_into_rigid_fragments()
    {
        var system = Create();
        int count = system.SpawnRigidBodyFromMesh(Glass(1.0f), 0.002f, breakSpeed: 2.0f, out int glass, mass: 0.001f);
        Assert.True(count > 1000);

        int shatterEvents = 0;
        system.GroupShattered += (_, _, _) => shatterEvents++;

        for (int i = 0; i < 90; i++) system.Update(Dt);

        Assert.False(system.IsGroupAlive(glass));
        Assert.True(shatterEvents >= 1);
        Assert.True(system.RigidGroupCount >= 10, $"only {system.RigidGroupCount} fragments");

        // Every particle is on or above the floor, and the pieces have spread out.
        float spread = 0;
        for (int i = 0; i < count; i++)
        {
            var p = system.GetParticle(i).Position;
            Assert.False(float.IsNaN(p.X));
            Assert.True(p.Y >= 0.0019f, $"particle below the floor at y={p.Y}");
            spread = MathF.Max(spread, MathF.Sqrt(p.X * p.X + p.Z * p.Z));
        }
        Assert.True(spread > 0.08f, "the fragments did not scatter");

        // A fragment stays rigid.
        int fragment = Enumerable.Range(0, count).Select(system.GetParticleGroup).First(g => g >= 0);
        var members = Enumerable.Range(0, count).Where(i => system.GetParticleGroup(i) == fragment).ToArray();
        var before = members.Select(i => system.GetParticle(i).Position).ToArray();
        for (int i = 0; i < 30; i++) system.Update(Dt);
        Assert.True(MaxPairDistanceError(system, members, before) < 1e-3f, "a fragment deformed");
    }

    [Fact]
    public void Glass_put_down_gently_does_not_break()
    {
        var system = Create();
        system.SpawnRigidBodyFromMesh(Glass(0.05f), 0.002f, breakSpeed: 2.0f, out int glass, mass: 0.001f);

        for (int i = 0; i < 120; i++) system.Update(Dt);

        Assert.True(system.IsGroupAlive(glass));
        Assert.Equal(1, system.RigidGroupCount);
        var (center, _) = system.GetGroupPose(glass);
        Assert.InRange(center.Y, 0.02f, 0.08f);
    }

    [Fact]
    public void Rigid_body_can_be_shattered_on_demand_and_compacted()
    {
        var system = Create();
        var box = TriangleMesh.CreateBox(new Vector3(-0.1f, 0.5f, -0.1f), new Vector3(0.1f, 0.7f, 0.1f));
        int count = system.SpawnRigidBodyFromMesh(box, 0.01f, float.PositiveInfinity, out int group);

        int fragments = system.ShatterGroup(group, new Vector3(0, 0.5f, 0));

        Assert.True(fragments >= 2);
        Assert.Equal(fragments, system.RigidGroupCount);

        system.Compact();
        Assert.Equal(count, system.Count);
        for (int i = 0; i < 60; i++) system.Update(Dt);
        Assert.Equal(fragments, system.RigidGroupCount);
    }
}
