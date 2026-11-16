using Artemis.Bodies;
using Artemis.Core;
using Artemis.Materials;
using Artemis.Simulation;

namespace Artemis.Tests;

public class SolverStressTests
{
    private const double Dt = 1.0 / 60.0;

    private static PhysicsMaterial Mat(double restitution, double friction)
        => new() { Restitution = restitution, StaticFriction = friction, DynamicFriction = friction };

    private static RigidBody AddIncline(PhysicsWorld world, double angle, double friction)
    {
        var ramp = RigidBody.CreateStaticBox(Vector3D.Zero, new Vector3D(20, 0.5, 5), Mat(0, friction));
        ramp.Rotation = Quaternion.FromAxisAngle(Vector3D.Forward, angle);
        world.AddBody(ramp);
        return ramp;
    }

    private static RigidBody AddBoxOnIncline(PhysicsWorld world, double angle, double friction)
    {
        var rot = Quaternion.FromAxisAngle(Vector3D.Forward, angle);
        var box = RigidBody.CreateBox(rot * new Vector3D(0, 1.0, 0), new Vector3D(0.5, 0.5, 0.5), 1.0, Mat(0, friction));
        box.Rotation = rot;
        box.LinearDamping = 0;
        world.AddBody(box);
        return box;
    }

    [Fact]
    public void Box_stays_put_on_an_incline_shallower_than_the_friction_angle()
    {
        // tan(20 deg) = 0.36 < mu = 0.6
        double angle = 20 * PhysicsConstants.DegToRad;
        var world = new PhysicsWorld();
        AddIncline(world, angle, 0.6);
        var box = AddBoxOnIncline(world, angle, 0.6);
        var start = box.Position;

        for (int i = 0; i < 240; i++) world.Step(Dt);

        Assert.True((box.Position - start).Magnitude < 0.02, $"box crept {(box.Position - start).Magnitude} m");
    }

    [Fact]
    public void Box_slides_down_a_steep_incline_with_the_analytic_acceleration()
    {
        // a = g (sin(theta) - mu cos(theta))
        const double mu = 0.2;
        double angle = 35 * PhysicsConstants.DegToRad;
        var world = new PhysicsWorld { AllowSleeping = false };
        AddIncline(world, angle, mu);
        var box = AddBoxOnIncline(world, angle, mu);
        var start = box.Position;

        const double time = 1.5;
        for (int i = 0; i < (int)(time / Dt); i++) world.Step(Dt);

        double a = PhysicsConstants.EarthGravity * (Math.Sin(angle) - mu * Math.Cos(angle));
        double expected = 0.5 * a * time * time;
        double travelled = (box.Position - start).Magnitude;
        Assert.InRange(travelled, expected * 0.9, expected * 1.1);
    }

    [Fact]
    public void Pile_of_spheres_in_a_container_settles_without_gaining_energy()
    {
        var world = new PhysicsWorld();
        world.CreateGround(20, 20);
        world.CreateWalls(6, 10, 6);

        var rng = new Random(7);
        var balls = new List<RigidBody>();
        for (int i = 0; i < 300; i++)
        {
            var ball = RigidBody.CreateSphere(
                new Vector3D(rng.NextDouble() * 4 - 2, 1 + i * 0.12, rng.NextDouble() * 4 - 2),
                0.25 + rng.NextDouble() * 0.1, 1.0, Mat(0.2, 0.5));
            balls.Add(ball);
            world.AddBody(ball);
        }

        double maxSpeed = 0;
        for (int i = 0; i < 900; i++)
        {
            world.Step(Dt);
            foreach (var b in balls) maxSpeed = Math.Max(maxSpeed, b.Velocity.Magnitude);
        }

        // Free fall from the top of the column is the fastest anything may legitimately get.
        double freeFallLimit = Math.Sqrt(2 * PhysicsConstants.EarthGravity * (1 + 300 * 0.12)) * 1.05;
        Assert.True(maxSpeed < freeFallLimit, $"a sphere reached {maxSpeed} m/s");

        foreach (var b in balls)
        {
            Assert.False(double.IsNaN(b.Position.X));
            Assert.InRange(b.Position.Y, 0.15, 10);
            Assert.InRange(b.Position.X, -3.2, 3.2);
            Assert.InRange(b.Position.Z, -3.2, 3.2);
        }

        double meanSpeed = balls.Average(b => b.Velocity.Magnitude);
        Assert.True(meanSpeed < 0.1, $"pile is still moving at {meanSpeed} m/s on average");
        Assert.True(balls.Count(b => b.IsSleeping) > balls.Count / 2, "most of a settled pile should be asleep");
    }

    [Fact]
    public void Simulation_is_deterministic_including_the_multi_threaded_path()
    {
        static List<Vector3D> Simulate()
        {
            var world = new PhysicsWorld();
            world.CreateGround(40, 40);
            var rng = new Random(123);
            var bodies = new List<RigidBody>();
            for (int i = 0; i < 400; i++)
            {
                var pos = new Vector3D(rng.NextDouble() * 6 - 3, 1 + i * 0.1, rng.NextDouble() * 6 - 3);
                RigidBody body = i % 2 == 0
                    ? RigidBody.CreateSphere(pos, 0.3, 1.0)
                    : RigidBody.CreateBox(pos, new Vector3D(0.3, 0.2, 0.25), 1.5);
                bodies.Add(body);
                world.AddBody(body);
            }

            for (int i = 0; i < 240; i++) world.Step(Dt);
            return bodies.Select(b => b.Position).ToList();
        }

        var first = Simulate();
        var second = Simulate();

        for (int i = 0; i < first.Count; i++)
        {
            Assert.Equal(first[i].X, second[i].X);
            Assert.Equal(first[i].Y, second[i].Y);
            Assert.Equal(first[i].Z, second[i].Z);
        }
    }

    [Fact]
    public void Box_pyramid_stays_standing()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var boxes = new List<(RigidBody Body, Vector3D Start)>();
        const int rows = 5;
        for (int row = 0; row < rows; row++)
        {
            for (int k = 0; k < rows - row; k++)
            {
                var pos = new Vector3D(k * 1.05 + row * 0.525, 0.5 + row * 1.0, 0);
                var box = RigidBody.CreateBox(pos, new Vector3D(0.5, 0.5, 0.5), 1.0, Mat(0, 0.6));
                boxes.Add((box, pos));
                world.AddBody(box);
            }
        }

        for (int i = 0; i < 600; i++) world.Step(Dt);

        foreach (var (body, start) in boxes)
            Assert.True((body.Position - start).Magnitude < 0.1,
                $"box moved {(body.Position - start).Magnitude} m from {start}");
    }

    [Fact]
    public void Heavy_box_on_light_box_does_not_sink_or_explode()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var light = RigidBody.CreateBox(new Vector3D(0, 0.5, 0), new Vector3D(0.5, 0.5, 0.5), 1.0, Mat(0, 0.6));
        var heavy = RigidBody.CreateBox(new Vector3D(0, 1.5, 0), new Vector3D(0.5, 0.5, 0.5), 20.0, Mat(0, 0.6));
        world.AddBody(light);
        world.AddBody(heavy);

        for (int i = 0; i < 600; i++) world.Step(Dt);

        Assert.InRange(light.Position.Y, 0.42, 0.55);
        Assert.InRange(heavy.Position.Y, 1.35, 1.55);
    }

    [Fact]
    public void Free_spinning_asymmetric_body_keeps_bounded_kinetic_energy()
    {
        var box = RigidBody.CreateBox(Vector3D.Zero, new Vector3D(0.2, 0.5, 1.0), 3.0);
        box.AngularDamping = 0;
        box.CanSleep = false;
        // Spin close to the unstable intermediate axis.
        box.AngularVelocity = new Vector3D(0.1, 6.0, 0.1);

        double Energy()
        {
            var r = Matrix3x3.FromQuaternion(box.Rotation);
            var local = r.Transposed * box.AngularVelocity;
            return 0.5 * Vector3D.Dot(local, box.InertiaTensor * local);
        }

        double e0 = Energy();
        double max = e0, min = e0;
        for (int i = 0; i < 3600; i++)
        {
            box.Integrate(Dt);
            double e = Energy();
            max = Math.Max(max, e);
            min = Math.Min(min, e);
        }

        Assert.InRange(max / e0, 0.9, 1.15);
        Assert.InRange(min / e0, 0.85, 1.1);
    }

    [Fact]
    public void Rolling_sphere_reaches_pure_rolling_under_friction()
    {
        var world = new PhysicsWorld { AllowSleeping = false };
        var ground = world.CreateGround();
        ground.Material = Mat(0, 0.5);
        var ball = RigidBody.CreateSphere(new Vector3D(0, 0.5, 0), 0.5, 1.0, Mat(0, 0.5));
        ball.LinearDamping = 0;
        ball.AngularDamping = 0;
        ball.Velocity = new Vector3D(7, 0, 0);
        world.AddBody(ball);

        for (int i = 0; i < 180; i++) world.Step(Dt);

        // A sliding solid sphere ends up rolling at 5/7 of its initial speed.
        Assert.InRange(ball.Velocity.X, 4.8, 5.2);
        // Rolling without slipping: v = w x r  ->  w_z = -v / r
        Assert.InRange(-ball.AngularVelocity.Z * 0.5, ball.Velocity.X * 0.97, ball.Velocity.X * 1.03);
    }
}
