using Artemis.Bodies;
using Artemis.Collision;
using Artemis.Core;
using Artemis.Materials;
using Artemis.Simulation;

namespace Artemis.Tests;

public class RigidBodyWorldTests
{
    private const double Dt = 1.0 / 60.0;

    private static void Run(PhysicsWorld world, double seconds)
    {
        int steps = (int)Math.Round(seconds / Dt);
        for (int i = 0; i < steps; i++)
            world.Step(Dt);
    }

    private static PhysicsMaterial Mat(double restitution, double friction)
        => new() { Restitution = restitution, StaticFriction = friction, DynamicFriction = friction };

    [Fact]
    public void Sphere_dropped_on_ground_comes_to_rest_on_the_surface()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var ball = RigidBody.CreateSphere(new Vector3D(0, 3, 0), 0.5, 1.0);
        world.AddBody(ball);

        Run(world, 5);

        Assert.InRange(ball.Position.Y, 0.48, 0.52);
        Assert.True(ball.Velocity.Magnitude < 0.05, $"still moving: {ball.Velocity}");
    }

    [Fact]
    public void Sphere_lands_on_ground_far_from_the_ground_centre()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var ball = RigidBody.CreateSphere(new Vector3D(23.3, 3, -17.1), 0.5, 1.0);
        world.AddBody(ball);

        Run(world, 5);

        Assert.InRange(ball.Position.Y, 0.48, 0.52);
    }

    [Fact]
    public void Box_resting_off_centre_on_ground_does_not_spin_or_drift()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var box = RigidBody.CreateBox(new Vector3D(10, 0.6, 5), new Vector3D(0.5, 0.5, 0.5), 2.0);
        world.AddBody(box);

        Run(world, 4);

        Assert.InRange(box.Position.Y, 0.48, 0.52);
        Assert.InRange(box.Position.X, 9.95, 10.05);
        Assert.InRange(box.Position.Z, 4.95, 5.05);
        Assert.True(Quaternion.Angle(box.Rotation, Quaternion.Identity) < 0.02,
            $"box rotated by {Quaternion.Angle(box.Rotation, Quaternion.Identity)} rad");
    }

    [Fact]
    public void Stack_of_boxes_stays_standing()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var boxes = new List<RigidBody>();
        for (int i = 0; i < 4; i++)
        {
            var box = RigidBody.CreateBox(new Vector3D(3, 0.5 + i * 1.0, 3), new Vector3D(0.5, 0.5, 0.5), 1.0,
                Mat(0, 0.6));
            boxes.Add(box);
            world.AddBody(box);
        }

        Run(world, 6);

        for (int i = 0; i < boxes.Count; i++)
        {
            Assert.InRange(boxes[i].Position.Y, 0.45 + i, 0.55 + i);
            Assert.InRange(boxes[i].Position.X, 2.9, 3.1);
            Assert.InRange(boxes[i].Position.Z, 2.9, 3.1);
        }
    }

    [Fact]
    public void Tilted_box_dropped_on_ground_settles_flat_on_a_face()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var box = RigidBody.CreateBox(new Vector3D(0, 2, 0), new Vector3D(0.5, 0.5, 0.5), 1.0, Mat(0.1, 0.6));
        box.Rotation = Quaternion.FromAxisAngle(new Vector3D(1, 0, 1), 0.5);
        world.AddBody(box);

        Run(world, 8);

        Assert.InRange(box.Position.Y, 0.47, 0.53);
        Assert.True(box.Velocity.Magnitude < 0.05, $"still moving: {box.Velocity}");
        // One local axis must be aligned with world up.
        double best = 0;
        foreach (var axis in new[] { Vector3D.Right, Vector3D.Up, Vector3D.Forward })
            best = Math.Max(best, Math.Abs(Vector3D.Dot(box.Rotation * axis, Vector3D.Up)));
        Assert.True(best > 0.999, $"box is not lying flat (alignment {best})");
    }

    [Fact]
    public void Sliding_box_is_stopped_by_coulomb_friction_over_the_expected_distance()
    {
        const double mu = 0.5;
        const double v0 = 4.0;
        var world = new PhysicsWorld();
        var ground = world.CreateGround();
        ground.Material = Mat(0, mu);
        var box = RigidBody.CreateBox(new Vector3D(0, 0.5, 0), new Vector3D(0.5, 0.25, 0.5), 1.0, Mat(0, mu));
        box.Position = new Vector3D(0, 0.25, 0);
        box.LinearDamping = 0;
        box.Velocity = new Vector3D(v0, 0, 0);
        world.AddBody(box);

        Run(world, 4);

        double expected = v0 * v0 / (2 * mu * PhysicsConstants.EarthGravity); // 1.63 m
        Assert.True(box.Velocity.Magnitude < 0.05, $"still sliding: {box.Velocity}");
        Assert.InRange(box.Position.X, expected * 0.85, expected * 1.15);
    }

    [Fact]
    public void Box_on_a_frictional_ground_does_not_slide_without_a_push()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var box = RigidBody.CreateBox(new Vector3D(-7, 0.5, 2), new Vector3D(0.5, 0.5, 0.5), 1.0);
        world.AddBody(box);

        Run(world, 3);

        Assert.InRange(box.Position.X, -7.01, -6.99);
        Assert.InRange(box.Position.Z, 1.99, 2.01);
    }

    [Fact]
    public void Resting_body_falls_asleep_and_is_woken_by_an_impact()
    {
        var world = new PhysicsWorld();
        world.CreateGround();
        var box = RigidBody.CreateBox(new Vector3D(0, 0.5, 0), new Vector3D(0.5, 0.5, 0.5), 1.0);
        world.AddBody(box);

        Run(world, 3);
        Assert.True(box.IsSleeping, "resting box never went to sleep");
        Assert.InRange(box.Position.Y, 0.48, 0.52);

        var ball = RigidBody.CreateSphere(new Vector3D(-3, 0.5, 0), 0.3, 1.0);
        ball.Velocity = new Vector3D(10, 0, 0);
        world.AddBody(ball);
        Run(world, 0.5);

        Assert.False(box.IsSleeping && box.Position.X < 0.01, "box was not woken by the impact");
        Assert.True(box.Position.X > 0.01, "box did not react to the impact");
    }

    [Fact]
    public void Head_on_elastic_collision_of_equal_spheres_swaps_velocities()
    {
        var world = new PhysicsWorld(Vector3D.Zero);
        var a = RigidBody.CreateSphere(new Vector3D(-2, 0, 0), 0.5, 1.0, Mat(1, 0));
        var b = RigidBody.CreateSphere(new Vector3D(2, 0, 0), 0.5, 1.0, Mat(1, 0));
        a.LinearDamping = b.LinearDamping = 0;
        a.Velocity = new Vector3D(3, 0, 0);
        world.AddBody(a);
        world.AddBody(b);

        Run(world, 2);

        Assert.InRange(a.Velocity.X, -0.05, 0.05);
        Assert.InRange(b.Velocity.X, 2.95, 3.05);
    }

    [Fact]
    public void Momentum_is_conserved_in_an_oblique_collision_of_unequal_spheres()
    {
        var world = new PhysicsWorld(Vector3D.Zero);
        var a = RigidBody.CreateSphere(new Vector3D(-2, 0.3, 0), 0.5, 1.0, Mat(0.5, 0.4));
        var b = RigidBody.CreateSphere(new Vector3D(2, 0, 0.2), 0.8, 5.0, Mat(0.5, 0.4));
        a.LinearDamping = b.LinearDamping = 0;
        a.Velocity = new Vector3D(6, 0, 0);
        b.Velocity = new Vector3D(-1, 0, 0);
        var before = a.Velocity * a.Mass + b.Velocity * b.Mass;
        world.AddBody(a);
        world.AddBody(b);

        Run(world, 2);

        var after = a.Velocity * a.Mass + b.Velocity * b.Mass;
        Assert.True((after - before).Magnitude < 1e-6, $"momentum changed by {(after - before).Magnitude}");
        Assert.True(a.Velocity.X < 5.9, "spheres never collided");
    }

    [Fact]
    public void Spheres_in_diagonally_adjacent_broadphase_cells_still_collide()
    {
        // Centres sit in cells (0,1,0) and (1,0,0) of a 2 m grid.
        var world = new PhysicsWorld(Vector3D.Zero);
        var a = RigidBody.CreateSphere(new Vector3D(1.7, 2.3, 0.5), 0.6, 1.0);
        var b = RigidBody.CreateSphere(new Vector3D(2.3, 1.7, 0.5), 0.6, 1.0);
        world.AddBody(a);
        world.AddBody(b);

        world.Step(Dt);

        Assert.Single(world.LastCollisions);
    }

    [Fact]
    public void Each_overlapping_pair_is_reported_exactly_once()
    {
        var world = new PhysicsWorld(Vector3D.Zero);
        // A large sphere straddling many cells against a small one.
        world.AddBody(RigidBody.CreateSphere(new Vector3D(0, 0, 0), 5, 1.0));
        world.AddBody(RigidBody.CreateSphere(new Vector3D(5.2, 0, 0), 0.5, 1.0));

        world.Step(Dt);

        Assert.Single(world.LastCollisions);
    }

    [Fact]
    public void Fast_sphere_does_not_tunnel_through_a_thin_wall_with_ccd()
    {
        var world = new PhysicsWorld(Vector3D.Zero) { UseCCD = true };
        var wall = RigidBody.CreateStaticBox(new Vector3D(5, 0, 0), new Vector3D(0.05, 5, 5));
        var bullet = RigidBody.CreateSphere(new Vector3D(0, 0, 0), 0.1, 0.05);
        bullet.LinearDamping = 0;
        bullet.Velocity = new Vector3D(300, 0, 0); // 5 m per step, wall is 0.1 m thick
        world.AddBody(wall);
        world.AddBody(bullet);

        Run(world, 1);

        Assert.True(bullet.Position.X < 5, $"bullet tunnelled to x={bullet.Position.X}");
    }

    [Fact]
    public void Sphere_bounces_to_roughly_restitution_squared_of_the_drop_height()
    {
        var world = new PhysicsWorld();
        var ground = world.CreateGround();
        ground.Material = Mat(0.8, 0.5);
        var ball = RigidBody.CreateSphere(new Vector3D(4, 2.5, 4), 0.5, 1.0, Mat(0.8, 0.5));
        ball.LinearDamping = 0;
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

        // Drop height 2.0 m, e = 0.8 -> rebound ~1.28 m.
        Assert.InRange(apex, 1.1, 1.4);
    }

    [Fact]
    public void Kinematic_body_moves_at_its_velocity_and_ignores_gravity_and_impacts()
    {
        var world = new PhysicsWorld();
        var platform = RigidBody.CreateBox(new Vector3D(0, 1, 0), new Vector3D(2, 0.1, 2), 1.0);
        platform.BodyType = BodyType.Kinematic;
        platform.Velocity = new Vector3D(1, 0, 0);
        var ball = RigidBody.CreateSphere(new Vector3D(0, 3, 0), 0.3, 10.0);
        world.AddBody(platform);
        world.AddBody(ball);

        Run(world, 1);

        Assert.InRange(platform.Position.X, 0.99, 1.01);
        Assert.InRange(platform.Position.Y, 0.999, 1.001);
        Assert.True(ball.Position.Y > 1.3, "ball fell through the kinematic platform");
    }
}

public class RigidBodyPropertyTests
{
    [Fact]
    public void Box_inertia_matches_the_solid_cuboid_formula()
    {
        var box = RigidBody.CreateBox(Vector3D.Zero, new Vector3D(1, 2, 3), 6.0);
        // I = m/3 * (b^2 + c^2) with half extents
        Assert.Equal(6.0 / 3 * (4 + 9), box.InertiaTensor.M00, 9);
        Assert.Equal(6.0 / 3 * (1 + 9), box.InertiaTensor.M11, 9);
        Assert.Equal(6.0 / 3 * (1 + 4), box.InertiaTensor.M22, 9);
    }

    [Fact]
    public void Sphere_inertia_uses_the_actual_radius()
    {
        var sphere = RigidBody.CreateSphere(Vector3D.Zero, 2.0, 5.0);
        Assert.Equal(0.4 * 5.0 * 4.0, sphere.InertiaTensor.M00, 9);
    }

    [Fact]
    public void Changing_shape_parameters_updates_inertia_and_bounds()
    {
        var body = new RigidBody(Vector3D.Zero, 2.0) { ShapeType = CollisionShapeType.Sphere };
        body.Radius = 3.0;
        Assert.Equal(0.4 * 2.0 * 9.0, body.InertiaTensor.M00, 9);
        Assert.Equal(3.0, body.BoundingBox.Max.X, 9);
    }

    [Fact]
    public void Static_and_kinematic_bodies_have_no_inverse_mass_or_inverse_inertia()
    {
        var ground = RigidBody.CreateStaticBox(Vector3D.Zero, new Vector3D(5, 1, 5));
        Assert.Equal(0, ground.InverseMass);
        Assert.Equal(Matrix3x3.Zero, ground.WorldInverseInertiaTensor);

        var body = RigidBody.CreateBox(Vector3D.Zero, Vector3D.One, 1.0);
        body.BodyType = BodyType.Kinematic;
        Assert.Equal(0, body.InverseMass);
        Assert.Equal(Matrix3x3.Zero, body.WorldInverseInertiaTensor);

        body.BodyType = BodyType.Dynamic;
        Assert.Equal(1.0, body.InverseMass, 9);
        Assert.NotEqual(Matrix3x3.Zero, body.WorldInverseInertiaTensor);
    }

    [Fact]
    public void Torque_free_spinning_body_keeps_its_angular_momentum()
    {
        var box = RigidBody.CreateBox(Vector3D.Zero, new Vector3D(0.2, 0.5, 1.0), 3.0);
        box.AngularDamping = 0;
        box.CanSleep = false;
        box.AngularVelocity = new Vector3D(1.0, 4.0, 0.3);

        Vector3D L(RigidBody b)
        {
            var r = Matrix3x3.FromQuaternion(b.Rotation);
            return r * (b.InertiaTensor * (r.Transposed * b.AngularVelocity));
        }

        var before = L(box);
        for (int i = 0; i < 600; i++)
            box.Integrate(1.0 / 60.0);
        var after = L(box);

        Assert.True((after - before).Magnitude / before.Magnitude < 0.02,
            $"angular momentum drifted from {before} to {after}");
    }
}

public class CollisionDetectorTests
{
    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Sphere_box_normal_points_from_body_a_to_body_b(bool sphereFirst)
    {
        var box = RigidBody.CreateStaticBox(new Vector3D(0, -0.5, 0), new Vector3D(5, 0.5, 5));
        var sphere = RigidBody.CreateSphere(new Vector3D(1, 0.4, 1), 0.5, 1.0);

        var info = sphereFirst ? CollisionDetector.Detect(sphere, box) : CollisionDetector.Detect(box, sphere);

        Assert.True(info.HasCollision);
        var aToB = info.BodyB.Position - info.BodyA.Position;
        Assert.True(Vector3D.Dot(info.Normal, new Vector3D(0, aToB.Y, 0)) > 0, $"normal {info.Normal} points B->A");
        Assert.Equal(0.1, info.Penetration, 6);
        Assert.Equal(0.0, info.ContactPoint.Y, 6);
    }

    [Fact]
    public void Sphere_centred_inside_box_still_gets_a_unit_normal()
    {
        var box = RigidBody.CreateStaticBox(Vector3D.Zero, Vector3D.One);
        var sphere = RigidBody.CreateSphere(Vector3D.Zero, 0.5, 1.0);

        var info = CollisionDetector.Detect(sphere, box);

        Assert.True(info.HasCollision);
        Assert.Equal(1.0, info.Normal.Magnitude, 6);
    }

    [Fact]
    public void Box_box_contact_lies_on_the_touching_faces_not_between_centres()
    {
        var ground = RigidBody.CreateStaticBox(new Vector3D(0, -0.5, 0), new Vector3D(50, 0.5, 50));
        var box = RigidBody.CreateBox(new Vector3D(20, 0.45, 10), new Vector3D(0.5, 0.5, 0.5), 1.0);

        var info = CollisionDetector.Detect(ground, box);

        Assert.True(info.HasCollision);
        Assert.Equal(0.05, info.Penetration, 6);
        Assert.True(info.Normal.Y > 0.999);
        Assert.InRange(info.ContactPoint.X, 19.5, 20.5);
        Assert.InRange(info.ContactPoint.Z, 9.5, 10.5);
        Assert.InRange(info.ContactPoint.Y, -0.06, 0.01);
    }

    [Fact]
    public void Box_box_face_contact_produces_a_four_point_manifold()
    {
        var ground = RigidBody.CreateStaticBox(new Vector3D(0, -0.5, 0), new Vector3D(50, 0.5, 50));
        var box = RigidBody.CreateBox(new Vector3D(20, 0.45, 10), new Vector3D(0.5, 0.5, 0.5), 1.0);

        var info = CollisionDetector.Detect(ground, box);

        Assert.Equal(4, info.ContactCount);
        for (int i = 0; i < info.ContactCount; i++)
        {
            var c = info.GetContact(i);
            Assert.Equal(0.5, Math.Abs(c.Position.X - 20), 6);
            Assert.Equal(0.5, Math.Abs(c.Position.Z - 10), 6);
            Assert.Equal(0.05, c.Penetration, 6);
        }
    }

    [Fact]
    public void Raycast_hits_the_nearest_body_and_reports_the_surface_normal()
    {
        var world = new PhysicsWorld(Vector3D.Zero);
        world.AddBody(RigidBody.CreateSphere(new Vector3D(5, 0, 0), 1, 1));
        world.AddBody(RigidBody.CreateBox(new Vector3D(10, 0, 0), Vector3D.One, 1));

        var hit = world.Raycast(Vector3D.Zero, new Vector3D(1, 0, 0));

        Assert.NotNull(hit);
        Assert.Equal(4.0, hit!.Value.Distance, 6);
        Assert.Equal(-1.0, hit.Value.Normal.X, 6);
    }
}
