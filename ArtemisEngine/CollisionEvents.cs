namespace ArtemisEngine;

public class CollisionEventArgs : EventArgs
{
    public RigidBody BodyA { get; set; }
    public RigidBody BodyB { get; set; }
    public Collision CollisionInfo { get; set; }

    public CollisionEventArgs(RigidBody bodyA, RigidBody bodyB, Collision collision)
    {
        BodyA = bodyA;
        BodyB = bodyB;
        CollisionInfo = collision;
    }
}

public delegate void CollisionEventHandler(object sender, CollisionEventArgs e);

public interface ICollisionListener
{
    void OnCollisionEnter(RigidBody other, Collision collision);
    void OnCollisionStay(RigidBody other, Collision collision);
    void OnCollisionExit(RigidBody other);
}
