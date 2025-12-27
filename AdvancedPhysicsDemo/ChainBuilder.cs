using ArtemisEngine;

namespace AdvancedPhysicsDemo;

/// <summary>
/// Creates a chain or rope using distance joints
/// Demonstrates joint constraints
/// </summary>
public class ChainBuilder
{
    public static List<RigidBody> CreateChain(PhysicsWorld world, Vector2 startPos, int segments, float segmentLength, bool isRope = false)
    {
        var bodies = new List<RigidBody>();
        var joints = new List<Joint>();

        // Create first segment (fixed)
        var firstShape = isRope ? (Shape)new CircleShape(0.3f) : new BoxShape(0.4f, segmentLength * 0.8f);
        var firstBody = new RigidBody(startPos, 0, firstShape, isStatic: true);
        firstBody.CollisionLayer = CollisionLayers.Static;
        bodies.Add(firstBody);
        world.AddBody(firstBody);

        RigidBody previousBody = firstBody;

        // Create chain segments
        for (int i = 1; i < segments; i++)
        {
            Vector2 position = startPos + new Vector2(0, -i * segmentLength);

            Shape shape = isRope ? (Shape)new CircleShape(0.3f) : new BoxShape(0.4f, segmentLength * 0.8f);
            float mass = isRope ? 1.0f : 2.0f;

            var body = new RigidBody(position, mass, shape);
            body.Friction = 0.5f;
            body.Restitution = 0.1f;
            body.CollisionLayer = CollisionLayers.Default;

            bodies.Add(body);
            world.AddBody(body);

            // Connect with previous segment
            var joint = new DistanceJoint(
                previousBody,
                body,
                new Vector2(0, isRope ? 0 : -segmentLength * 0.4f),
                new Vector2(0, isRope ? 0 : segmentLength * 0.4f)
            );
            joint.Stiffness = isRope ? 0.8f : 0.95f;

            joints.Add(joint);
            world.AddJoint(joint);

            previousBody = body;
        }

        return bodies;
    }

    public static List<RigidBody> CreateBridge(PhysicsWorld world, Vector2 startPos, int segments, float segmentWidth)
    {
        var bodies = new List<RigidBody>();

        // Create fixed anchor on left
        var leftAnchor = new RigidBody(startPos, 0, new BoxShape(0.5f, 0.5f), isStatic: true);
        leftAnchor.CollisionLayer = CollisionLayers.Static;
        bodies.Add(leftAnchor);
        world.AddBody(leftAnchor);

        RigidBody previousBody = leftAnchor;

        // Create bridge segments
        for (int i = 1; i <= segments; i++)
        {
            Vector2 position = startPos + new Vector2(i * segmentWidth, 0);

            var shape = new BoxShape(segmentWidth * 0.9f, 0.3f);
            var body = new RigidBody(position, 3.0f, shape);
            body.Friction = 0.6f;
            body.Restitution = 0.1f;

            bodies.Add(body);
            world.AddBody(body);

            // Connect with rope joints for flexibility
            var joint = new RopeJoint(previousBody, body, segmentWidth * 1.1f);
            joint.Stiffness = 0.9f;
            world.AddJoint(joint);

            previousBody = body;
        }

        // Create fixed anchor on right
        Vector2 rightPos = startPos + new Vector2((segments + 1) * segmentWidth, 0);
        var rightAnchor = new RigidBody(rightPos, 0, new BoxShape(0.5f, 0.5f), isStatic: true);
        rightAnchor.CollisionLayer = CollisionLayers.Static;
        bodies.Add(rightAnchor);
        world.AddBody(rightAnchor);

        var finalJoint = new RopeJoint(previousBody, rightAnchor, segmentWidth * 1.1f);
        finalJoint.Stiffness = 0.9f;
        world.AddJoint(finalJoint);

        return bodies;
    }
}
