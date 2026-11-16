using Artemis.Physics2D;

namespace Artemis.Tests;

public class World2DTests
{
    private const double Dt = 1.0 / 60.0;

    private static void Run(PhysicsWorld2D world, double seconds)
    {
        for (int i = 0; i < (int)Math.Round(seconds / Dt); i++)
            world.Step(Dt);
    }

    private static PhysicsWorld2D WorldWithGround(out RigidBody2D ground)
    {
        var world = new PhysicsWorld2D();
        ground = RigidBody2D.CreateStaticBox(new Vector2D(0, -0.5), 100, 1);
        world.AddBody(ground);
        return world;
    }

    [Fact]
    public void Circle_dropped_on_ground_rests_on_the_surface()
    {
        var world = WorldWithGround(out _);
        var ball = RigidBody2D.CreateCircle(new Vector2D(13, 3), 0.5);
        world.AddBody(ball);

        Run(world, 5);

        Assert.InRange(ball.Position.Y, 0.47, 0.52);
        Assert.True(ball.Velocity.Magnitude < 0.05, $"still moving {ball.Velocity}");
    }

    [Fact]
    public void Box_dropped_on_ground_rests_flat_without_spinning()
    {
        var world = WorldWithGround(out _);
        var box = RigidBody2D.CreateBox(new Vector2D(-21, 2), 1, 1);
        world.AddBody(box);

        Run(world, 5);

        Assert.InRange(box.Position.Y, 0.47, 0.52);
        Assert.InRange(box.Position.X, -21.05, -20.95);
        Assert.InRange(box.Rotation, -0.02, 0.02);
    }

    [Fact]
    public void Free_fall_matches_the_analytic_solution()
    {
        var world = new PhysicsWorld2D(new Vector2D(0, -10));
        var ball = RigidBody2D.CreateCircle(new Vector2D(0, 100), 0.5);
        world.AddBody(ball);

        Run(world, 2);

        Assert.InRange(100 - ball.Position.Y, 19.8, 20.2);
        Assert.Equal(-20.0, ball.Velocity.Y, 6);
    }

    [Fact]
    public void Stack_of_boxes_stays_standing()
    {
        var world = WorldWithGround(out _);
        var boxes = new List<RigidBody2D>();
        for (int i = 0; i < 6; i++)
        {
            var box = RigidBody2D.CreateBox(new Vector2D(2, 0.5 + i), 1, 1);
            box.Friction = 0.6;
            box.Restitution = 0;
            boxes.Add(box);
            world.AddBody(box);
        }

        Run(world, 8);

        for (int i = 0; i < boxes.Count; i++)
        {
            Assert.InRange(boxes[i].Position.Y, 0.4 + i, 0.55 + i);
            Assert.InRange(boxes[i].Position.X, 1.9, 2.1);
        }
    }

    [Fact]
    public void Sliding_box_stops_after_the_coulomb_distance()
    {
        const double mu = 0.5, v0 = 4.0;
        var world = WorldWithGround(out var ground);
        ground.Friction = mu;
        var box = RigidBody2D.CreateBox(new Vector2D(0, 0.25), 1, 0.5);
        box.Friction = mu;
        box.Restitution = 0;
        box.Velocity = new Vector2D(v0, 0);
        world.AddBody(box);

        Run(world, 4);

        double expected = v0 * v0 / (2 * mu * 9.81);
        Assert.True(box.Velocity.Magnitude < 0.05, $"still sliding {box.Velocity}");
        Assert.InRange(box.Position.X, expected * 0.85, expected * 1.15);
    }

    [Fact]
    public void Head_on_elastic_collision_of_equal_circles_swaps_velocities()
    {
        var world = new PhysicsWorld2D(Vector2D.Zero);
        var a = RigidBody2D.CreateCircle(new Vector2D(-2, 0), 0.5);
        var b = RigidBody2D.CreateCircle(new Vector2D(2, 0), 0.5);
        a.Restitution = b.Restitution = 1;
        a.Friction = b.Friction = 0;
        a.Velocity = new Vector2D(3, 0);
        world.AddBody(a);
        world.AddBody(b);

        Run(world, 2);

        Assert.InRange(a.Velocity.X, -0.05, 0.05);
        Assert.InRange(b.Velocity.X, 2.95, 3.05);
    }

    [Fact]
    public void Momentum_is_conserved_in_an_oblique_collision()
    {
        var world = new PhysicsWorld2D(Vector2D.Zero);
        var a = RigidBody2D.CreateCircle(new Vector2D(-2, 0.3), 0.5, 1.0);
        var b = RigidBody2D.CreateCircle(new Vector2D(2, 0), 0.8, 3.0);
        a.Velocity = new Vector2D(6, 0);
        b.Velocity = new Vector2D(-1, 0);
        var before = a.Velocity * a.Mass + b.Velocity * b.Mass;
        world.AddBody(a);
        world.AddBody(b);

        Run(world, 2);

        var after = a.Velocity * a.Mass + b.Velocity * b.Mass;
        Assert.True((after - before).Magnitude < 1e-6 * before.Magnitude + 1e-9, $"momentum changed by {(after - before).Magnitude}");
        Assert.True(a.Velocity.X < 5.9, "circles never collided");
    }

    [Fact]
    public void Ball_bounces_to_restitution_squared_of_the_drop_height()
    {
        var world = WorldWithGround(out var ground);
        ground.Restitution = 0.8;
        var ball = RigidBody2D.CreateCircle(new Vector2D(4, 2.5), 0.5);
        ball.Restitution = 0.8;
        world.AddBody(ball);

        bool hit = false;
        double apex = 0;
        for (int i = 0; i < 180; i++)
        {
            world.Step(Dt);
            if (!hit && ball.Velocity.Y > 0) hit = true;
            if (hit)
            {
                apex = Math.Max(apex, ball.Position.Y - 0.5);
                if (ball.Velocity.Y < 0) break;
            }
        }

        Assert.InRange(apex, 1.1, 1.4);
    }

    [Fact]
    public void Resting_box_falls_asleep_and_wakes_on_impact()
    {
        var world = WorldWithGround(out _);
        var box = RigidBody2D.CreateBox(new Vector2D(0, 0.5), 1, 1);
        world.AddBody(box);

        Run(world, 3);
        Assert.True(box.IsSleeping, "resting box never went to sleep");

        var ball = RigidBody2D.CreateCircle(new Vector2D(-3, 0.5), 0.3);
        ball.Velocity = new Vector2D(10, 0);
        world.AddBody(ball);
        Run(world, 0.5);

        Assert.True(box.Position.X > 0.01, "box did not react to the impact");
    }

    [Fact]
    public void Fast_bullet_does_not_tunnel_through_a_thin_wall()
    {
        var world = new PhysicsWorld2D(Vector2D.Zero);
        var wall = RigidBody2D.CreateStaticBox(new Vector2D(5, 0), 0.1, 10);
        var bullet = RigidBody2D.CreateCircle(new Vector2D(0.3, 0), 0.1);
        bullet.IsBullet = true;
        bullet.Velocity = new Vector2D(300, 0);
        world.AddBody(wall);
        world.AddBody(bullet);

        Run(world, 1);

        Assert.True(bullet.Position.X < 5, $"bullet tunnelled to x={bullet.Position.X}");
    }

    [Fact]
    public void Pendulum_on_a_distance_joint_keeps_its_length_and_never_gains_energy()
    {
        const double g = 9.81;
        var world = new PhysicsWorld2D(new Vector2D(0, -g));
        var anchor = RigidBody2D.CreateStaticBox(new Vector2D(0, 5), 0.2, 0.2);
        var bob = RigidBody2D.CreateCircle(new Vector2D(2, 5), 0.2);
        bob.CanSleep = false;
        world.AddBody(anchor);
        world.AddBody(bob);
        world.CreateDistanceJoint(anchor, bob, anchor.Position, bob.Position);

        double Energy() => g * bob.Position.Y + 0.5 * bob.Velocity.MagnitudeSquared;
        double initial = Energy();

        double maxLengthError = 0, lowest = double.MaxValue, maxEnergy = initial;
        for (int i = 0; i < 600; i++)
        {
            world.Step(Dt);
            maxLengthError = Math.Max(maxLengthError, Math.Abs((bob.Position - anchor.Position).Magnitude - 2));
            lowest = Math.Min(lowest, bob.Position.Y);
            maxEnergy = Math.Max(maxEnergy, Energy());
        }

        Assert.True(maxLengthError < 0.01, $"rod stretched by {maxLengthError}");
        Assert.InRange(lowest, 2.99, 3.01);
        Assert.True(maxEnergy <= initial * 1.001, $"energy grew from {initial} to {maxEnergy}");
        // The first-order integrator is slightly dissipative, but the swing must not die.
        Assert.True(Energy() > g * 4.0, "pendulum lost most of its energy");
    }

    [Theory]
    [InlineData("polygon")]
    [InlineData("box")]
    [InlineData("circle")]
    public void Bodies_rest_on_edge_and_chain_ground(string kind)
    {
        var world = new PhysicsWorld2D();
        world.AddBody(RigidBody2D.CreateEdge(new Vector2D(-10, 0), new Vector2D(0, 0)));
        world.AddBody(RigidBody2D.CreateChain(new[] { new Vector2D(0, 0), new Vector2D(5, 0), new Vector2D(10, 0) }));

        RigidBody2D Make(Vector2D pos) => kind switch
        {
            "polygon" => new RigidBody2D(pos, 1.0, PolygonShape.CreateRegular(6, 0.5)),
            "box" => RigidBody2D.CreateBox(pos, 1, 1),
            _ => RigidBody2D.CreateCircle(pos, 0.5)
        };

        var onEdge = Make(new Vector2D(-5, 2));
        var onChain = Make(new Vector2D(6, 2));
        world.AddBody(onEdge);
        world.AddBody(onChain);

        Run(world, 5);

        // A regular hexagon of circumradius 0.5 lies on a face at height 0.433.
        double lo = kind == "polygon" ? 0.40 : 0.47;
        Assert.InRange(onEdge.Position.Y, lo, 0.52);
        Assert.InRange(onChain.Position.Y, lo, 0.52);
    }

    [Fact]
    public void Polygon_box_and_circle_collide_with_each_other()
    {
        var world = WorldWithGround(out _);
        var box = RigidBody2D.CreateBox(new Vector2D(0, 0.5), 2, 1);
        var hexagon = new RigidBody2D(new Vector2D(0, 2), 1.0, PolygonShape.CreateRegular(6, 0.5));
        var ball = RigidBody2D.CreateCircle(new Vector2D(0, 3.5), 0.3);
        world.AddBody(box);
        world.AddBody(hexagon);
        world.AddBody(ball);

        Run(world, 6);

        Assert.InRange(box.Position.Y, 0.47, 0.52);
        Assert.True(hexagon.Position.Y > 1.35, $"hexagon sank into the box (y={hexagon.Position.Y})");
        Assert.True(ball.Position.Y > 0.28, "ball fell through");
        Assert.True(ball.Position.Y > 1.2 || Math.Abs(ball.Position.X) > 0.9, "ball passed through the hexagon or box");
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Circle_box_normal_points_from_body_a_to_body_b(bool circleFirst)
    {
        var box = RigidBody2D.CreateStaticBox(new Vector2D(0, -0.5), 10, 1);
        var circle = RigidBody2D.CreateCircle(new Vector2D(1, 0.4), 0.5);
        var manifold = new Manifold2D();

        bool hit = circleFirst
            ? CollisionDetector2D.DetectCollision(circle, box, manifold)
            : CollisionDetector2D.DetectCollision(box, circle, manifold);

        Assert.True(hit);
        var aToB = manifold.BodyB.Position - manifold.BodyA.Position;
        Assert.True(manifold.Contacts[0].Normal.Y * aToB.Y > 0, "normal points from B to A");
        Assert.Equal(0.1, manifold.Contacts[0].Depth, 6);
    }

    [Fact]
    public void Box_resting_on_box_gets_a_two_point_manifold()
    {
        var ground = RigidBody2D.CreateStaticBox(new Vector2D(0, -0.5), 10, 1);
        var box = RigidBody2D.CreateBox(new Vector2D(3, 0.45), 1, 1);
        var manifold = new Manifold2D();

        Assert.True(CollisionDetector2D.DetectCollision(ground, box, manifold));
        Assert.Equal(2, manifold.ContactCount);
        Assert.Equal(0.05, manifold.Contacts[0].Depth, 6);
        Assert.Equal(0.05, manifold.Contacts[1].Depth, 6);
        Assert.Equal(1.0, Math.Abs(manifold.Contacts[0].Point.X - manifold.Contacts[1].Point.X), 6);
    }

    [Fact]
    public void Body_made_dynamic_after_creation_as_static_can_move()
    {
        var world = new PhysicsWorld2D();
        var body = RigidBody2D.CreateStaticBox(new Vector2D(0, 10), 1, 1);
        body.Mass = 2;
        body.Inertia = 1;
        world.AddBody(body);

        body.BodyType = BodyType2D.Dynamic;
        Run(world, 1);

        Assert.True(body.InverseMass > 0);
        Assert.True(body.Position.Y < 6, "body never started to fall");
    }

    [Fact]
    public void Simulation_is_deterministic()
    {
        static List<Vector2D> Simulate()
        {
            var world = new PhysicsWorld2D();
            world.AddBody(RigidBody2D.CreateStaticBox(new Vector2D(0, -0.5), 100, 1));
            var rng = new Random(5);
            var bodies = new List<RigidBody2D>();
            for (int i = 0; i < 150; i++)
            {
                var pos = new Vector2D(rng.NextDouble() * 6 - 3, 1 + i * 0.4);
                var body = i % 2 == 0 ? RigidBody2D.CreateCircle(pos, 0.3) : RigidBody2D.CreateBox(pos, 0.6, 0.4);
                bodies.Add(body);
                world.AddBody(body);
            }
            for (int i = 0; i < 300; i++) world.Step(Dt);
            return bodies.Select(b => b.Position).ToList();
        }

        var first = Simulate();
        var second = Simulate();
        for (int i = 0; i < first.Count; i++)
        {
            Assert.Equal(first[i].X, second[i].X);
            Assert.Equal(first[i].Y, second[i].Y);
        }
    }

    [Fact]
    public void Pile_of_bodies_settles_inside_its_container()
    {
        var world = WorldWithGround(out _);
        world.AddBody(RigidBody2D.CreateStaticBox(new Vector2D(-4.5, 10), 1, 20));
        world.AddBody(RigidBody2D.CreateStaticBox(new Vector2D(4.5, 10), 1, 20));
        var rng = new Random(11);
        var bodies = new List<RigidBody2D>();
        for (int i = 0; i < 120; i++)
        {
            var pos = new Vector2D(rng.NextDouble() * 6 - 3, 1 + i * 0.5);
            var body = i % 2 == 0 ? RigidBody2D.CreateCircle(pos, 0.3) : RigidBody2D.CreateBox(pos, 0.6, 0.5);
            bodies.Add(body);
            world.AddBody(body);
        }

        Run(world, 20);

        foreach (var b in bodies)
        {
            Assert.False(double.IsNaN(b.Position.X));
            Assert.InRange(b.Position.X, -4.1, 4.1);
            Assert.InRange(b.Position.Y, 0.15, 30);
        }
        Assert.True(bodies.Average(b => b.Velocity.Magnitude) < 0.1, "pile is still moving");
    }
}

public class Joints2DTests
{
    private const double Dt = 1.0 / 60.0;

    [Fact]
    public void Hanging_chain_of_revolute_joints_stays_connected_and_comes_to_rest()
    {
        var world = new PhysicsWorld2D();
        var anchor = RigidBody2D.CreateStaticBox(new Vector2D(0, 10), 0.2, 0.2);
        world.AddBody(anchor);

        var links = new List<RigidBody2D>();
        var previous = anchor;
        for (int i = 0; i < 12; i++)
        {
            // Laid out horizontally so the chain has to swing down.
            var link = RigidBody2D.CreateBox(new Vector2D(0.25 + i * 0.5, 10), 0.5, 0.1);
            link.LinearDamping = 0.5;
            link.AngularDamping = 0.5;
            world.AddBody(link);
            world.CreateRevoluteJoint(previous, link, new Vector2D(i * 0.5, 10));
            links.Add(link);
            previous = link;
        }

        double maxGap = 0;
        for (int step = 0; step < 1800; step++)
        {
            world.Step(Dt);
            for (int i = 1; i < links.Count; i++)
            {
                var endOfPrevious = links[i - 1].LocalToWorld(new Vector2D(0.25, 0));
                var startOfThis = links[i].LocalToWorld(new Vector2D(-0.25, 0));
                maxGap = Math.Max(maxGap, (endOfPrevious - startOfThis).Magnitude);
            }
        }

        Assert.True(maxGap < 0.1, $"chain links separated by {maxGap} m");
        var last = links[^1];
        Assert.False(double.IsNaN(last.Position.X));
        // 12 links of 0.5 m hanging straight down from y = 10.
        Assert.InRange(last.Position.Y, 4.0, 4.6);
        Assert.InRange(last.Position.X, -0.3, 0.3);
    }

    [Fact]
    public void Projectile_knocks_over_a_tower_and_everything_settles()
    {
        var world = new PhysicsWorld2D();
        world.AddBody(RigidBody2D.CreateStaticBox(new Vector2D(0, -0.5), 200, 1));

        var tower = new List<RigidBody2D>();
        for (int i = 0; i < 8; i++)
        {
            var block = RigidBody2D.CreateBox(new Vector2D(10, 0.5 + i), 1, 1);
            tower.Add(block);
            world.AddBody(block);
        }

        for (int i = 0; i < 120; i++) world.Step(Dt);
        Assert.All(tower, b => Assert.InRange(b.Position.X, 9.95, 10.05));
        Assert.All(tower, b => Assert.True(b.IsSleeping));

        var shot = RigidBody2D.CreateCircle(new Vector2D(0, 5.5), 0.5, 2.0);
        shot.IsBullet = true;
        shot.Velocity = new Vector2D(15, 2);
        world.AddBody(shot);

        for (int i = 0; i < 1200; i++) world.Step(Dt);

        Assert.True(tower.Count(b => Math.Abs(b.Position.X - 10) > 0.5) >= 3, "tower was not knocked over");
        foreach (var b in tower)
        {
            Assert.False(double.IsNaN(b.Position.X));
            Assert.True(b.Position.Y > 0.4, $"block sank into the ground (y={b.Position.Y})");
            Assert.True(b.Velocity.Magnitude < 0.2, $"block still moving at {b.Velocity.Magnitude}");
        }
    }
}
