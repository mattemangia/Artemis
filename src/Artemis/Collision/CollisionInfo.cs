using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Collision
{
    /// <summary>
    /// Contains information about a collision between two bodies.
    /// </summary>
    public struct CollisionInfo
    {
        /// <summary>
        /// First body involved in the collision.
        /// </summary>
        public IPhysicsBody BodyA;

        /// <summary>
        /// Second body involved in the collision.
        /// </summary>
        public IPhysicsBody BodyB;

        /// <summary>
        /// Collision normal (pointing from A to B).
        /// </summary>
        public Vector3D Normal;

        /// <summary>
        /// Penetration depth (positive if overlapping).
        /// </summary>
        public double Penetration;

        /// <summary>
        /// Contact point in world coordinates.
        /// </summary>
        public Vector3D ContactPoint;

        /// <summary>
        /// Whether a collision was detected.
        /// </summary>
        public bool HasCollision;

        /// <summary>
        /// Creates a collision info indicating no collision.
        /// </summary>
        public static CollisionInfo NoCollision => new() { HasCollision = false };

        /// <summary>
        /// Creates a collision info with the specified parameters.
        /// </summary>
        public static CollisionInfo Create(
            IPhysicsBody bodyA,
            IPhysicsBody bodyB,
            Vector3D normal,
            double penetration,
            Vector3D contactPoint)
        {
            return new CollisionInfo
            {
                BodyA = bodyA,
                BodyB = bodyB,
                Normal = normal,
                Penetration = penetration,
                ContactPoint = contactPoint,
                HasCollision = true
            };
        }

        /// <summary>
        /// Returns a collision info with bodies swapped.
        /// </summary>
        public readonly CollisionInfo Swapped()
        {
            return new CollisionInfo
            {
                BodyA = BodyB,
                BodyB = BodyA,
                Normal = -Normal,
                Penetration = Penetration,
                ContactPoint = ContactPoint,
                HasCollision = HasCollision
            };
        }

        public override readonly string ToString()
            => HasCollision
                ? $"Collision(Normal: {Normal}, Depth: {Penetration:F4}, Contact: {ContactPoint})"
                : "NoCollision";
    }

    /// <summary>
    /// Contains contact manifold data for persistent contact tracking.
    /// </summary>
    public class ContactManifold
    {
        /// <summary>
        /// Maximum number of contact points in a manifold.
        /// </summary>
        public const int MaxContacts = 4;

        /// <summary>
        /// First body.
        /// </summary>
        public IPhysicsBody BodyA { get; set; } = null!;

        /// <summary>
        /// Second body.
        /// </summary>
        public IPhysicsBody BodyB { get; set; } = null!;

        /// <summary>
        /// Contact points.
        /// </summary>
        public ContactPoint[] Contacts { get; } = new ContactPoint[MaxContacts];

        /// <summary>
        /// Number of active contact points.
        /// </summary>
        public int ContactCount { get; set; }

        /// <summary>
        /// Average normal of all contacts.
        /// </summary>
        public Vector3D Normal { get; set; }

        /// <summary>
        /// Adds a contact point to the manifold.
        /// </summary>
        public void AddContact(ContactPoint contact)
        {
            if (ContactCount < MaxContacts)
            {
                Contacts[ContactCount++] = contact;
            }
            else
            {
                // Replace the contact with the smallest penetration
                int minIndex = 0;
                double minPen = Contacts[0].Penetration;
                for (int i = 1; i < MaxContacts; i++)
                {
                    if (Contacts[i].Penetration < minPen)
                    {
                        minPen = Contacts[i].Penetration;
                        minIndex = i;
                    }
                }

                if (contact.Penetration > minPen)
                {
                    Contacts[minIndex] = contact;
                }
            }
        }

        /// <summary>
        /// Clears all contacts.
        /// </summary>
        public void Clear()
        {
            ContactCount = 0;
        }
    }

    /// <summary>
    /// Represents a single contact point.
    /// </summary>
    public struct ContactPoint
    {
        /// <summary>
        /// Position in world space.
        /// </summary>
        public Vector3D Position;

        /// <summary>
        /// Position relative to body A.
        /// </summary>
        public Vector3D LocalA;

        /// <summary>
        /// Position relative to body B.
        /// </summary>
        public Vector3D LocalB;

        /// <summary>
        /// Contact normal.
        /// </summary>
        public Vector3D Normal;

        /// <summary>
        /// Penetration depth.
        /// </summary>
        public double Penetration;

        /// <summary>
        /// Accumulated normal impulse (for warm starting).
        /// </summary>
        public double NormalImpulse;

        /// <summary>
        /// Accumulated tangent impulse (friction).
        /// </summary>
        public double TangentImpulse;

        /// <summary>
        /// Lifetime of this contact in frames.
        /// </summary>
        public int Lifetime;
    }
}
