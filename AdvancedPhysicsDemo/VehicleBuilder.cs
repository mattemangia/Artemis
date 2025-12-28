
namespace AdvancedPhysicsDemo;

/// <summary>
/// Creates a simple vehicle with wheels connected via revolute joints
/// Demonstrates revolute joints and spring suspension
/// </summary>
public class VehicleBuilder
{
    public class Vehicle
    {
        public RigidBody Body { get; set; }
        public RigidBody FrontWheel { get; set; }
        public RigidBody BackWheel { get; set; }
        public List<Joint> Joints { get; set; } = new();

        public Vehicle(RigidBody body, RigidBody frontWheel, RigidBody backWheel)
        {
            Body = body;
            FrontWheel = frontWheel;
            BackWheel = backWheel;
        }

        public void ApplyAcceleration(float force)
        {
            FrontWheel.ApplyTorque(force);
            BackWheel.ApplyTorque(force);
        }
    }

    public static Vehicle CreateVehicle(PhysicsWorld world, Vector2 position)
    {
        // Create car body
        var bodyShape = new BoxShape(2.0f, 0.75f);
        var carBody = new RigidBody(position, 20.0f, bodyShape);
        carBody.Friction = 0.5f;
        carBody.Restitution = 0.2f;
        carBody.CollisionLayer = CollisionLayers.Player;
        world.AddBody(carBody);

        // Create front wheel
        var wheelShape = new CircleShape(0.8f);
        var frontWheel = new RigidBody(position + new Vector2(1.5f, -1.2f), 5.0f, wheelShape);
        frontWheel.Friction = 1.0f;
        frontWheel.Restitution = 0.3f;
        frontWheel.CollisionLayer = CollisionLayers.Player;
        world.AddBody(frontWheel);

        // Create back wheel
        var backWheel = new RigidBody(position + new Vector2(-1.5f, -1.2f), 5.0f, wheelShape);
        backWheel.Friction = 1.0f;
        backWheel.Restitution = 0.3f;
        backWheel.CollisionLayer = CollisionLayers.Player;
        world.AddBody(backWheel);

        // Connect wheels with spring joints for suspension
        var frontSpring = new SpringJoint(carBody, frontWheel, 1.2f, 200f, 10f);
        frontSpring.LocalAnchorA = new Vector2(1.5f, -0.5f);
        frontSpring.LocalAnchorB = new Vector2(0, 0);
        world.AddJoint(frontSpring);

        var backSpring = new SpringJoint(carBody, backWheel, 1.2f, 200f, 10f);
        backSpring.LocalAnchorA = new Vector2(-1.5f, -0.5f);
        backSpring.LocalAnchorB = new Vector2(0, 0);
        world.AddJoint(backSpring);

        var vehicle = new Vehicle(carBody, frontWheel, backWheel);
        vehicle.Joints.Add(frontSpring);
        vehicle.Joints.Add(backSpring);

        return vehicle;
    }
}
