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
        [ThreadStatic]
        private static ContactSolver? _solver;

        [ThreadStatic]
        private static CollisionInfo[]? _single;

        /// <summary>
        /// Resolves a single collision between two bodies: applies the normal and
        /// friction impulses and pushes the bodies out of overlap.
        /// <para>
        /// This looks at one pair in isolation. To resolve all collisions of a step
        /// together (required for stable stacking) use <see cref="ContactSolver"/>,
        /// which is what the physics worlds do.
        /// </para>
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

            // An explicit resolve is an impact: sleeping bodies take part in it.
            if (bodyA.BodyType == BodyType.Dynamic) bodyA.WakeUp();
            if (bodyB.BodyType == BodyType.Dynamic) bodyB.WakeUp();

            var solver = _solver ??= new ContactSolver { PositionCorrectionFactor = 0.4, WarmStarting = false };
            var single = _single ??= new CollisionInfo[1];

            solver.RestitutionThreshold = restitutionThreshold;
            single[0] = collision;
            solver.SolveVelocities(single);
            solver.SolvePositions();
            single[0] = default;
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
