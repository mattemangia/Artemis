using System;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Collision event arguments for 2D collisions.
    /// </summary>
    public class CollisionEvent2D : EventArgs
    {
        public RigidBody2D BodyA { get; }
        public RigidBody2D BodyB { get; }
        public Manifold2D Manifold { get; }

        public CollisionEvent2D(RigidBody2D bodyA, RigidBody2D bodyB, Manifold2D manifold)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            Manifold = manifold;
        }
    }
}
