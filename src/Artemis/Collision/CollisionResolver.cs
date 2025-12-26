using System;
using Artemis.Bodies;
using Artemis.Core;
using Artemis.Materials;

namespace Artemis.Collision
{
    /// <summary>
    /// Resolves collisions between physics bodies using impulse-based methods.
    /// </summary>
    public static class CollisionResolver
    {
        /// <summary>
        /// Resolves a collision between two bodies.
        /// </summary>
        /// <param name="collision">The collision information.</param>
        /// <param name="restitutionThreshold">Velocity threshold below which restitution is zero.</param>
        public static void Resolve(CollisionInfo collision, double restitutionThreshold = 1.0)
        {
            if (!collision.HasCollision)
                return;

            var bodyA = collision.BodyA;
            var bodyB = collision.BodyB;

            // Skip if both bodies are static/kinematic
            if (bodyA.BodyType != BodyType.Dynamic && bodyB.BodyType != BodyType.Dynamic)
                return;

            // Separate the bodies to resolve penetration
            SeparateBodies(collision);

            // Apply collision impulse
            ApplyCollisionImpulse(collision, restitutionThreshold);

            // Apply friction impulse
            ApplyFrictionImpulse(collision);

            // Wake up sleeping bodies
            bodyA.WakeUp();
            bodyB.WakeUp();
        }

        /// <summary>
        /// Separates overlapping bodies.
        /// </summary>
        private static void SeparateBodies(CollisionInfo collision)
        {
            var bodyA = collision.BodyA;
            var bodyB = collision.BodyB;
            var normal = collision.Normal;
            double penetration = collision.Penetration;

            // Small bias to prevent jitter
            double slop = 0.001;
            double percent = 0.8; // Usually 20-80%

            double correction = Math.Max(penetration - slop, 0) * percent;

            double totalInverseMass = bodyA.InverseMass + bodyB.InverseMass;
            if (totalInverseMass <= PhysicsConstants.Epsilon)
                return;

            var correctionVector = normal * (correction / totalInverseMass);

            if (bodyA.BodyType == BodyType.Dynamic)
                bodyA.Position -= correctionVector * bodyA.InverseMass;

            if (bodyB.BodyType == BodyType.Dynamic)
                bodyB.Position += correctionVector * bodyB.InverseMass;
        }

        /// <summary>
        /// Applies the collision impulse for velocity resolution.
        /// </summary>
        private static void ApplyCollisionImpulse(CollisionInfo collision, double restitutionThreshold)
        {
            var bodyA = collision.BodyA;
            var bodyB = collision.BodyB;
            var normal = collision.Normal;
            var contactPoint = collision.ContactPoint;

            // Get relative velocity at contact point
            var rA = contactPoint - bodyA.Position;
            var rB = contactPoint - bodyB.Position;

            var velA = bodyA.Velocity + Vector3D.Cross(bodyA.AngularVelocity, rA);
            var velB = bodyB.Velocity + Vector3D.Cross(bodyB.AngularVelocity, rB);
            var relativeVelocity = velB - velA;

            // Relative velocity along normal
            double velAlongNormal = Vector3D.Dot(relativeVelocity, normal);

            // Do not resolve if velocities are separating
            if (velAlongNormal > 0)
                return;

            // Calculate restitution (average of materials)
            double restitution = (bodyA.Material.Restitution + bodyB.Material.Restitution) * 0.5;

            // Disable restitution for slow collisions to prevent jitter
            if (Math.Abs(velAlongNormal) < restitutionThreshold)
                restitution = 0;

            // Calculate impulse scalar
            double j = -(1 + restitution) * velAlongNormal;
            j /= CalculateImpulseDenominator(bodyA, bodyB, rA, rB, normal);

            // Apply impulse
            var impulse = normal * j;

            if (bodyA.BodyType == BodyType.Dynamic)
            {
                bodyA.ApplyImpulse(-impulse);
                if (bodyA is RigidBody rbA)
                {
                    var angImpulse = Vector3D.Cross(rA, -impulse);
                    bodyA.AngularVelocity += rbA.WorldInverseInertiaTensor * angImpulse;
                }
            }

            if (bodyB.BodyType == BodyType.Dynamic)
            {
                bodyB.ApplyImpulse(impulse);
                if (bodyB is RigidBody rbB)
                {
                    var angImpulse = Vector3D.Cross(rB, impulse);
                    bodyB.AngularVelocity += rbB.WorldInverseInertiaTensor * angImpulse;
                }
            }
        }

        /// <summary>
        /// Applies friction impulse.
        /// </summary>
        private static void ApplyFrictionImpulse(CollisionInfo collision)
        {
            var bodyA = collision.BodyA;
            var bodyB = collision.BodyB;
            var normal = collision.Normal;
            var contactPoint = collision.ContactPoint;

            var rA = contactPoint - bodyA.Position;
            var rB = contactPoint - bodyB.Position;

            var velA = bodyA.Velocity + Vector3D.Cross(bodyA.AngularVelocity, rA);
            var velB = bodyB.Velocity + Vector3D.Cross(bodyB.AngularVelocity, rB);
            var relativeVelocity = velB - velA;

            // Get tangent velocity
            var tangent = relativeVelocity - normal * Vector3D.Dot(relativeVelocity, normal);
            double tangentSpeed = tangent.Magnitude;

            if (tangentSpeed < PhysicsConstants.Epsilon)
                return;

            tangent = tangent / tangentSpeed;

            // Calculate friction coefficients (geometric mean)
            double staticFriction = Math.Sqrt(bodyA.Material.StaticFriction * bodyB.Material.StaticFriction);
            double dynamicFriction = Math.Sqrt(bodyA.Material.DynamicFriction * bodyB.Material.DynamicFriction);

            // Calculate tangent impulse
            double jt = -tangentSpeed;
            jt /= CalculateImpulseDenominator(bodyA, bodyB, rA, rB, tangent);

            // Get normal impulse magnitude for Coulomb friction
            double velAlongNormal = Vector3D.Dot(relativeVelocity, normal);
            double jn = Math.Abs(velAlongNormal) /
                       CalculateImpulseDenominator(bodyA, bodyB, rA, rB, normal);

            // Coulomb's law: use static or dynamic friction
            Vector3D frictionImpulse;
            if (Math.Abs(jt) < jn * staticFriction)
            {
                // Static friction
                frictionImpulse = tangent * jt;
            }
            else
            {
                // Dynamic friction
                frictionImpulse = tangent * (-jn * dynamicFriction);
            }

            // Apply friction impulse
            if (bodyA.BodyType == BodyType.Dynamic)
            {
                bodyA.ApplyImpulse(-frictionImpulse);
                if (bodyA is RigidBody rbA)
                {
                    var angImpulse = Vector3D.Cross(rA, -frictionImpulse);
                    bodyA.AngularVelocity += rbA.WorldInverseInertiaTensor * angImpulse;
                }
            }

            if (bodyB.BodyType == BodyType.Dynamic)
            {
                bodyB.ApplyImpulse(frictionImpulse);
                if (bodyB is RigidBody rbB)
                {
                    var angImpulse = Vector3D.Cross(rB, frictionImpulse);
                    bodyB.AngularVelocity += rbB.WorldInverseInertiaTensor * angImpulse;
                }
            }
        }

        /// <summary>
        /// Calculates the denominator for impulse computation.
        /// </summary>
        private static double CalculateImpulseDenominator(
            IPhysicsBody bodyA,
            IPhysicsBody bodyB,
            Vector3D rA,
            Vector3D rB,
            Vector3D direction)
        {
            double denom = bodyA.InverseMass + bodyB.InverseMass;

            if (bodyA is RigidBody rbA)
            {
                var crossA = Vector3D.Cross(rA, direction);
                var rotA = rbA.WorldInverseInertiaTensor * crossA;
                denom += Vector3D.Dot(Vector3D.Cross(rotA, rA), direction);
            }

            if (bodyB is RigidBody rbB)
            {
                var crossB = Vector3D.Cross(rB, direction);
                var rotB = rbB.WorldInverseInertiaTensor * crossB;
                denom += Vector3D.Dot(Vector3D.Cross(rotB, rB), direction);
            }

            return Math.Max(denom, PhysicsConstants.Epsilon);
        }

        /// <summary>
        /// Resolves a collision with plasticity effects.
        /// </summary>
        /// <param name="collision">The collision information.</param>
        /// <param name="impactEnergy">Output: the energy of the impact.</param>
        public static void ResolveWithPlasticity(CollisionInfo collision, out double impactEnergy)
        {
            impactEnergy = 0;

            if (!collision.HasCollision)
                return;

            var bodyA = collision.BodyA;
            var bodyB = collision.BodyB;

            // Calculate relative velocity before resolution
            var relVelBefore = bodyB.Velocity - bodyA.Velocity;
            double speedBefore = relVelBefore.Magnitude;

            // Standard resolution
            Resolve(collision);

            // Calculate relative velocity after resolution
            var relVelAfter = bodyB.Velocity - bodyA.Velocity;
            double speedAfter = relVelAfter.Magnitude;

            // Calculate energy dissipation
            double reducedMass = 1.0 / (bodyA.InverseMass + bodyB.InverseMass);
            impactEnergy = 0.5 * reducedMass * (speedBefore * speedBefore - speedAfter * speedAfter);

            // Apply plastic deformation if materials support it
            if (impactEnergy > 0)
            {
                ApplyPlasticDeformation(bodyA, impactEnergy * 0.5);
                ApplyPlasticDeformation(bodyB, impactEnergy * 0.5);
            }
        }

        private static void ApplyPlasticDeformation(IPhysicsBody body, double energy)
        {
            if (body.Material.Ductility <= 0 || body.Material.YieldStrength <= 0)
                return;

            // Calculate stress from impact energy
            // This is a simplified model
            double volume = 1.0; // Would need actual volume calculation
            double stress = Math.Sqrt(2 * energy * body.Material.YoungsModulus / volume);

            if (stress > body.Material.YieldStrength)
            {
                // Body would experience plastic deformation
                // In a full implementation, this could modify the shape
                // For now, we just note that it happened
                double plasticStrain = (stress - body.Material.YieldStrength) / body.Material.YoungsModulus;

                // Could trigger an event or modify body properties
                // body.OnPlasticDeformation?.Invoke(plasticStrain);
            }
        }
    }
}
