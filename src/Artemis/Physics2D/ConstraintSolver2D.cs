using System;
using System.Collections.Generic;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Contact manifold - stores persistent contact information between frames.
    /// Used for warm starting to improve solver convergence.
    /// </summary>
    public class ContactManifold2D
    {
        public RigidBody2D BodyA { get; set; }
        public RigidBody2D BodyB { get; set; }
        public List<ContactPoint2D> ContactPoints { get; set; } = new();
        public int FrameCount { get; set; } = 0; // How many frames this contact has persisted

        public class ContactPoint2D
        {
            public Vector2D LocalPointA { get; set; } // Contact point in body A's local space
            public Vector2D LocalPointB { get; set; } // Contact point in body B's local space
            public double NormalImpulse { get; set; } // Accumulated normal impulse
            public double TangentImpulse { get; set; } // Accumulated friction impulse
            public double Separation { get; set; } // Penetration depth (negative if penetrating)
            public Vector2D Normal { get; set; }
        }

        public string GetKey()
        {
            // Create unique key for this contact pair (order-independent)
            int hashA = BodyA.GetHashCode();
            int hashB = BodyB.GetHashCode();
            return hashA < hashB ? $"{hashA}_{hashB}" : $"{hashB}_{hashA}";
        }
    }

    /// <summary>
    /// Contact persistence - maintains contact information across frames.
    /// Enables warm starting for faster convergence.
    /// </summary>
    public class ContactPersistence2D
    {
        private Dictionary<string, ContactManifold2D> _manifolds = new();
        private const int MaxFrameAge = 3;
        private const double ContactMatchThreshold = 0.1;

        public void BeginFrame()
        {
            var keysToRemove = new List<string>();
            foreach (var pair in _manifolds)
            {
                pair.Value.FrameCount++;
                if (pair.Value.FrameCount > MaxFrameAge)
                {
                    keysToRemove.Add(pair.Key);
                }
            }

            foreach (var key in keysToRemove)
            {
                _manifolds.Remove(key);
            }
        }

        public ContactManifold2D GetOrCreateManifold(RigidBody2D bodyA, RigidBody2D bodyB)
        {
            var tempManifold = new ContactManifold2D { BodyA = bodyA, BodyB = bodyB };
            string key = tempManifold.GetKey();

            if (!_manifolds.TryGetValue(key, out var manifold))
            {
                manifold = new ContactManifold2D { BodyA = bodyA, BodyB = bodyB };
                _manifolds[key] = manifold;
            }

            manifold.FrameCount = 0;
            return manifold;
        }

        public void UpdateManifold(ContactManifold2D manifold, Vector2D worldContactPoint, Vector2D normal, double separation)
        {
            Vector2D localA = WorldToLocal(manifold.BodyA, worldContactPoint);
            Vector2D localB = WorldToLocal(manifold.BodyB, worldContactPoint);

            ContactManifold2D.ContactPoint2D? matched = null;
            foreach (var cp in manifold.ContactPoints)
            {
                Vector2D worldA = LocalToWorld(manifold.BodyA, cp.LocalPointA);
                double dist = (worldA - worldContactPoint).Magnitude;

                if (dist < ContactMatchThreshold)
                {
                    matched = cp;
                    break;
                }
            }

            if (matched != null)
            {
                matched.LocalPointA = localA;
                matched.LocalPointB = localB;
                matched.Separation = separation;
                matched.Normal = normal;
            }
            else
            {
                manifold.ContactPoints.Add(new ContactManifold2D.ContactPoint2D
                {
                    LocalPointA = localA,
                    LocalPointB = localB,
                    NormalImpulse = 0,
                    TangentImpulse = 0,
                    Separation = separation,
                    Normal = normal
                });
            }
        }

        private Vector2D WorldToLocal(RigidBody2D body, Vector2D worldPoint)
        {
            Vector2D relative = worldPoint - body.Position;
            double cos = Math.Cos(-body.Rotation);
            double sin = Math.Sin(-body.Rotation);
            return new Vector2D(
                relative.X * cos - relative.Y * sin,
                relative.X * sin + relative.Y * cos
            );
        }

        private Vector2D LocalToWorld(RigidBody2D body, Vector2D localPoint)
        {
            double cos = Math.Cos(body.Rotation);
            double sin = Math.Sin(body.Rotation);
            return new Vector2D(
                localPoint.X * cos - localPoint.Y * sin + body.Position.X,
                localPoint.X * sin + localPoint.Y * cos + body.Position.Y
            );
        }

        public void Clear()
        {
            _manifolds.Clear();
        }
    }

    /// <summary>
    /// Iterative constraint solver using Sequential Impulse method.
    /// Much more stable than single-pass impulse resolution.
    /// </summary>
    public class SequentialImpulseSolver2D
    {
        public int VelocityIterations { get; set; } = 8;
        public int PositionIterations { get; set; } = 3;
        public double Baumgarte { get; set; } = 0.2; // Position correction factor
        public double SlowdownFactor { get; set; } = 1.0;

        private ContactPersistence2D _persistence = new();

        public void BeginFrame()
        {
            _persistence.BeginFrame();
        }

        /// <summary>
        /// Solve collision using iterative Sequential Impulse method with warm starting.
        /// </summary>
        public void SolveCollision(RigidBody2D bodyA, RigidBody2D bodyB, Manifold2D collision)
        {
            if (bodyA.BodyType == BodyType2D.Static && bodyB.BodyType == BodyType2D.Static)
                return;

            Vector2D normal = collision.ContactCount > 0 ? collision.Contacts[0].Normal : collision.Normal;
            Vector2D contactPoint = collision.ContactCount > 0
                ? collision.Contacts[0].Point
                : (bodyA.Position + bodyB.Position) * 0.5;
            double penetration = collision.Penetration;

            var manifold = _persistence.GetOrCreateManifold(bodyA, bodyB);
            _persistence.UpdateManifold(manifold, contactPoint, normal, -penetration);

            Vector2D rA = contactPoint - bodyA.Position;
            Vector2D rB = contactPoint - bodyB.Position;

            double rnA = Vector2D.Cross(rA, normal);
            double rnB = Vector2D.Cross(rB, normal);
            double effectiveMassNormal = bodyA.InverseMass + bodyB.InverseMass +
                                        rnA * rnA * bodyA.InverseInertia +
                                        rnB * rnB * bodyB.InverseInertia;

            if (effectiveMassNormal < 0.0001)
                return;

            effectiveMassNormal = 1.0 / effectiveMassNormal;

            Vector2D tangent = new Vector2D(-normal.Y, normal.X);

            double rtA = Vector2D.Cross(rA, tangent);
            double rtB = Vector2D.Cross(rB, tangent);
            double effectiveMassTangent = bodyA.InverseMass + bodyB.InverseMass +
                                         rtA * rtA * bodyA.InverseInertia +
                                         rtB * rtB * bodyB.InverseInertia;

            effectiveMassTangent = effectiveMassTangent < 0.0001 ? 0 : 1.0 / effectiveMassTangent;

            double restitution = Math.Max(bodyA.Restitution, bodyB.Restitution);
            double friction = Math.Sqrt(bodyA.Friction * bodyB.Friction);

            double accumulatedNormalImpulse = 0;
            double accumulatedTangentImpulse = 0;

            if (manifold.ContactPoints.Count > 0)
            {
                var cp = manifold.ContactPoints[0];
                accumulatedNormalImpulse = cp.NormalImpulse;
                accumulatedTangentImpulse = cp.TangentImpulse;
            }

            // Warm start
            if (accumulatedNormalImpulse > 0 || Math.Abs(accumulatedTangentImpulse) > 0)
            {
                Vector2D impulse = normal * accumulatedNormalImpulse + tangent * accumulatedTangentImpulse;

                if (bodyA.BodyType == BodyType2D.Dynamic)
                {
                    bodyA.Velocity -= impulse * bodyA.InverseMass;
                    bodyA.AngularVelocity -= Vector2D.Cross(rA, impulse) * bodyA.InverseInertia;
                }

                if (bodyB.BodyType == BodyType2D.Dynamic)
                {
                    bodyB.Velocity += impulse * bodyB.InverseMass;
                    bodyB.AngularVelocity += Vector2D.Cross(rB, impulse) * bodyB.InverseInertia;
                }
            }

            // Velocity iterations
            for (int i = 0; i < VelocityIterations; i++)
            {
                Vector2D vA = bodyA.Velocity + new Vector2D(-rA.Y * bodyA.AngularVelocity, rA.X * bodyA.AngularVelocity);
                Vector2D vB = bodyB.Velocity + new Vector2D(-rB.Y * bodyB.AngularVelocity, rB.X * bodyB.AngularVelocity);
                Vector2D relativeVelocity = vB - vA;

                double normalVelocity = Vector2D.Dot(relativeVelocity, normal);
                double restitutionBias = normalVelocity < -1.0 ? -restitution * normalVelocity : 0;

                double normalImpulseMagnitude = (-normalVelocity + restitutionBias) * effectiveMassNormal;

                double newAccumulatedNormal = Math.Max(accumulatedNormalImpulse + normalImpulseMagnitude, 0);
                normalImpulseMagnitude = newAccumulatedNormal - accumulatedNormalImpulse;
                accumulatedNormalImpulse = newAccumulatedNormal;

                Vector2D normalImpulse = normal * normalImpulseMagnitude;

                if (bodyA.BodyType == BodyType2D.Dynamic)
                {
                    bodyA.Velocity -= normalImpulse * bodyA.InverseMass;
                    bodyA.AngularVelocity -= Vector2D.Cross(rA, normalImpulse) * bodyA.InverseInertia;
                }

                if (bodyB.BodyType == BodyType2D.Dynamic)
                {
                    bodyB.Velocity += normalImpulse * bodyB.InverseMass;
                    bodyB.AngularVelocity += Vector2D.Cross(rB, normalImpulse) * bodyB.InverseInertia;
                }

                // Friction impulse
                if (effectiveMassTangent > 0 && friction > 0)
                {
                    vA = bodyA.Velocity + new Vector2D(-rA.Y * bodyA.AngularVelocity, rA.X * bodyA.AngularVelocity);
                    vB = bodyB.Velocity + new Vector2D(-rB.Y * bodyB.AngularVelocity, rB.X * bodyB.AngularVelocity);
                    relativeVelocity = vB - vA;

                    double tangentVelocity = Vector2D.Dot(relativeVelocity, tangent);
                    double tangentImpulseMagnitude = -tangentVelocity * effectiveMassTangent;

                    double maxFriction = friction * accumulatedNormalImpulse;
                    double newAccumulatedTangent = Math.Max(-maxFriction, Math.Min(accumulatedTangentImpulse + tangentImpulseMagnitude, maxFriction));
                    tangentImpulseMagnitude = newAccumulatedTangent - accumulatedTangentImpulse;
                    accumulatedTangentImpulse = newAccumulatedTangent;

                    Vector2D frictionImpulse = tangent * tangentImpulseMagnitude;

                    if (bodyA.BodyType == BodyType2D.Dynamic)
                    {
                        bodyA.Velocity -= frictionImpulse * bodyA.InverseMass;
                        bodyA.AngularVelocity -= Vector2D.Cross(rA, frictionImpulse) * bodyA.InverseInertia;
                    }

                    if (bodyB.BodyType == BodyType2D.Dynamic)
                    {
                        bodyB.Velocity += frictionImpulse * bodyB.InverseMass;
                        bodyB.AngularVelocity += Vector2D.Cross(rB, frictionImpulse) * bodyB.InverseInertia;
                    }
                }
            }

            // Store accumulated impulses
            if (manifold.ContactPoints.Count > 0)
            {
                manifold.ContactPoints[0].NormalImpulse = accumulatedNormalImpulse;
                manifold.ContactPoints[0].TangentImpulse = accumulatedTangentImpulse;
            }

            // Position iterations
            for (int i = 0; i < PositionIterations; i++)
            {
                if (penetration <= 0)
                    break;

                double correction = Math.Max(penetration - 0.01, 0) * Baumgarte;
                Vector2D positionImpulse = normal * correction;

                if (bodyA.BodyType == BodyType2D.Dynamic)
                {
                    double massRatioA = bodyA.InverseMass / (bodyA.InverseMass + bodyB.InverseMass);
                    bodyA.Position -= positionImpulse * massRatioA * SlowdownFactor;
                }

                if (bodyB.BodyType == BodyType2D.Dynamic)
                {
                    double massRatioB = bodyB.InverseMass / (bodyA.InverseMass + bodyB.InverseMass);
                    bodyB.Position += positionImpulse * massRatioB * SlowdownFactor;
                }

                penetration *= 0.5;
            }
        }

        public void Clear()
        {
            _persistence.Clear();
        }
    }
}
