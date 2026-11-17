using Artemis.Bodies;
using Artemis.Core;
using Artemis.Materials;
using Artemis.Simulation;

namespace Artemis.Tests;

public class OptimizedWorldTests
{
    private const double Dt = 1.0 / 60.0;

    private static PhysicsMaterial Mat(double restitution, double friction)
        => new() { Restitution = restitution, StaticFriction = friction, DynamicFriction = friction };

    private static void Run(OptimizedPhysicsWorld world, double seconds)
    {
        for (int i = 0; i < (int)Math.Round(seconds / Dt); i++)
            world.Step(Dt);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Bodies_fall_under_gravity_and_rest_on_the_ground(bool multiThreaded)
    {
        var world = new OptimizedPhysicsWorld { UseMultiThreading = multiThreaded };
        world.CreateGround();
        var ball = RigidBody.CreateSphere(new Vector3D(12, 3, -9), 0.5, 1.0);
        var box = RigidBody.CreateBox(new Vector3D(-15, 2, 20), new Vector3D(0.5, 0.5, 0.5), 2.0);
        world.AddBody(ball);
        world.AddBody(box);

        Run(world, 5);

        Assert.InRange(ball.Position.Y, 0.48, 0.52);
        Assert.InRange(box.Position.Y, 0.48, 0.52);
        Assert.True(Quaternion.Angle(box.Rotation, Quaternion.Identity) < 0.02);
        Assert.True(ball.IsSleeping && box.IsSleeping, "resting bodies should be asleep");
    }

    [Fact]
    public void Free_fall_matches_the_analytic_solution()
    {
        var world = new OptimizedPhysicsWorld(new Vector3D(0, -10, 0));
        var ball = RigidBody.CreateSphere(new Vector3D(0, 100, 0), 0.5, 1.0);
        ball.LinearDamping = 0;
        world.AddBody(ball);

        Run(world, 2);

        // Semi-implicit Euler overshoots 0.5 g t^2 by g t dt / 2.
        Assert.InRange(100 - ball.Position.Y, 20.0, 20.2);
        Assert.Equal(-20.0, ball.Velocity.Y, 6);
    }

    [Fact]
    public void Stack_of_boxes_stays_standing()
    {
        var world = new OptimizedPhysicsWorld();
        world.CreateGround();
        var boxes = new List<RigidBody>();
        for (int i = 0; i < 5; i++)
        {
            var box = RigidBody.CreateBox(new Vector3D(3, 0.5 + i, 3), new Vector3D(0.5, 0.5, 0.5), 1.0, Mat(0, 0.6));
            boxes.Add(box);
            world.AddBody(box);
        }

        Run(world, 6);

        for (int i = 0; i < boxes.Count; i++)
        {
            Assert.InRange(boxes[i].Position.Y, 0.45 + i, 0.55 + i);
            Assert.InRange(boxes[i].Position.X, 2.9, 3.1);
        }
    }

    [Fact]
    public void Sleeping_stack_is_woken_when_hit()
    {
        var world = new OptimizedPhysicsWorld();
        world.CreateGround();
        var box = RigidBody.CreateBox(new Vector3D(0, 0.5, 0), new Vector3D(0.5, 0.5, 0.5), 1.0);
        world.AddBody(box);
        Run(world, 3);
        Assert.True(box.IsSleeping);

        var ball = RigidBody.CreateSphere(new Vector3D(-3, 0.5, 0), 0.3, 1.0);
        ball.Velocity = new Vector3D(10, 0, 0);
        world.AddBody(ball);
        Run(world, 0.5);

        Assert.True(box.Position.X > 0.01, "box did not react to the impact");
    }

    [Fact]
    public void Removing_the_support_wakes_what_was_resting_on_it()
    {
        var world = new OptimizedPhysicsWorld();
        world.CreateGround();
        var lower = RigidBody.CreateBox(new Vector3D(0, 0.5, 0), new Vector3D(0.5, 0.5, 0.5), 1.0);
        var upper = RigidBody.CreateBox(new Vector3D(0, 1.5, 0), new Vector3D(0.5, 0.5, 0.5), 1.0);
        world.AddBody(lower);
        world.AddBody(upper);
        Run(world, 3);
        Assert.True(upper.IsSleeping);

        world.RemoveBody(lower);
        Run(world, 3);

        Assert.InRange(upper.Position.Y, 0.45, 0.55);
    }

    [Fact]
    public void Fast_sphere_does_not_tunnel_through_a_thin_wall()
    {
        var world = new OptimizedPhysicsWorld(Vector3D.Zero);
        var wall = RigidBody.CreateStaticBox(new Vector3D(5, 0, 0), new Vector3D(0.05, 5, 5));
        var bullet = RigidBody.CreateSphere(Vector3D.Zero, 0.1, 0.05);
        bullet.LinearDamping = 0;
        bullet.Velocity = new Vector3D(300, 0, 0);
        world.AddBody(wall);
        world.AddBody(bullet);

        Run(world, 1);

        Assert.True(bullet.Position.X < 5, $"bullet tunnelled to x={bullet.Position.X}");
    }

    [Fact]
    public void Body_switched_to_dynamic_after_being_added_starts_falling()
    {
        var world = new OptimizedPhysicsWorld();
        var body = RigidBody.CreateBox(new Vector3D(0, 10, 0), Vector3D.One, 1.0);
        body.BodyType = BodyType.Static;
        world.AddBody(body);
        Run(world, 0.5);
        Assert.Equal(10.0, body.Position.Y, 9);

        body.BodyType = BodyType.Dynamic;
        Run(world, 1);

        Assert.True(body.Position.Y < 6, "body never started to fall");
    }
}
