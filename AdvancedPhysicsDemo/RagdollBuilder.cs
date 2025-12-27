using ArtemisEngine;

namespace AdvancedPhysicsDemo;

/// <summary>
/// Creates a simple ragdoll using distance and revolute joints
/// Demonstrates character physics
/// </summary>
public class RagdollBuilder
{
    public class Ragdoll
    {
        public RigidBody Head { get; set; }
        public RigidBody Torso { get; set; }
        public RigidBody LeftArm { get; set; }
        public RigidBody RightArm { get; set; }
        public RigidBody LeftLeg { get; set; }
        public RigidBody RightLeg { get; set; }
        public List<Joint> Joints { get; set; } = new();

        public Ragdoll(RigidBody head, RigidBody torso, RigidBody leftArm, RigidBody rightArm,
                       RigidBody leftLeg, RigidBody rightLeg)
        {
            Head = head;
            Torso = torso;
            LeftArm = leftArm;
            RightArm = rightArm;
            LeftLeg = leftLeg;
            RightLeg = rightLeg;
        }

        public void ApplyImpulseToAll(Vector2 impulse)
        {
            Head.ApplyImpulse(impulse);
            Torso.ApplyImpulse(impulse);
            LeftArm.ApplyImpulse(impulse);
            RightArm.ApplyImpulse(impulse);
            LeftLeg.ApplyImpulse(impulse);
            RightLeg.ApplyImpulse(impulse);
        }
    }

    public static Ragdoll CreateRagdoll(PhysicsWorld world, Vector2 position)
    {
        // Head
        var head = new RigidBody(position + new Vector2(0, 2.5f), 2.0f, new CircleShape(0.4f));
        head.Friction = 0.5f;
        head.Restitution = 0.3f;
        head.CollisionLayer = CollisionLayers.Enemy;
        world.AddBody(head);

        // Torso
        var torso = new RigidBody(position + new Vector2(0, 1.2f), 5.0f, new BoxShape(1.0f, 1.5f));
        torso.Friction = 0.5f;
        torso.Restitution = 0.2f;
        torso.CollisionLayer = CollisionLayers.Enemy;
        world.AddBody(torso);

        // Left Arm
        var leftArm = new RigidBody(position + new Vector2(-0.8f, 1.5f), 1.5f, new BoxShape(0.3f, 1.2f));
        leftArm.Friction = 0.5f;
        leftArm.Restitution = 0.2f;
        leftArm.CollisionLayer = CollisionLayers.Enemy;
        world.AddBody(leftArm);

        // Right Arm
        var rightArm = new RigidBody(position + new Vector2(0.8f, 1.5f), 1.5f, new BoxShape(0.3f, 1.2f));
        rightArm.Friction = 0.5f;
        rightArm.Restitution = 0.2f;
        rightArm.CollisionLayer = CollisionLayers.Enemy;
        world.AddBody(rightArm);

        // Left Leg
        var leftLeg = new RigidBody(position + new Vector2(-0.3f, 0.0f), 2.0f, new BoxShape(0.4f, 1.5f));
        leftLeg.Friction = 0.6f;
        leftLeg.Restitution = 0.2f;
        leftLeg.CollisionLayer = CollisionLayers.Enemy;
        world.AddBody(leftLeg);

        // Right Leg
        var rightLeg = new RigidBody(position + new Vector2(0.3f, 0.0f), 2.0f, new BoxShape(0.4f, 1.5f));
        rightLeg.Friction = 0.6f;
        rightLeg.Restitution = 0.2f;
        rightLeg.CollisionLayer = CollisionLayers.Enemy;
        world.AddBody(rightLeg);

        var ragdoll = new Ragdoll(head, torso, leftArm, rightArm, leftLeg, rightLeg);

        // Connect head to torso (neck)
        var neck = new DistanceJoint(head, torso, new Vector2(0, -0.4f), new Vector2(0, 0.7f));
        neck.Stiffness = 0.8f;
        world.AddJoint(neck);
        ragdoll.Joints.Add(neck);

        // Connect left arm to torso (shoulder)
        var leftShoulder = new RevoluteJoint(leftArm, torso, position + new Vector2(-0.5f, 1.8f));
        leftShoulder.Stiffness = 0.7f;
        world.AddJoint(leftShoulder);
        ragdoll.Joints.Add(leftShoulder);

        // Connect right arm to torso (shoulder)
        var rightShoulder = new RevoluteJoint(rightArm, torso, position + new Vector2(0.5f, 1.8f));
        rightShoulder.Stiffness = 0.7f;
        world.AddJoint(rightShoulder);
        ragdoll.Joints.Add(rightShoulder);

        // Connect left leg to torso (hip)
        var leftHip = new RevoluteJoint(leftLeg, torso, position + new Vector2(-0.3f, 0.7f));
        leftHip.Stiffness = 0.8f;
        world.AddJoint(leftHip);
        ragdoll.Joints.Add(leftHip);

        // Connect right leg to torso (hip)
        var rightHip = new RevoluteJoint(rightLeg, torso, position + new Vector2(0.3f, 0.7f));
        rightHip.Stiffness = 0.8f;
        world.AddJoint(rightHip);
        ragdoll.Joints.Add(rightHip);

        return ragdoll;
    }
}
