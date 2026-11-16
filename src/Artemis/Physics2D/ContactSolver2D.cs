using System;
using System.Runtime.CompilerServices;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Sequential-impulse contact solver for <see cref="Manifold2D"/>.
    /// <para>
    /// Impulses are accumulated per contact point and clamped as a total, which is
    /// what makes resting contact, stacking and Coulomb friction behave. The
    /// accumulated impulses are kept in the contacts so the next step can start from
    /// them (warm starting).
    /// </para>
    /// <para>
    /// Per step: call <see cref="PreStep"/> once for every manifold after velocities
    /// have been integrated, then <see cref="SolveVelocity"/> for a number of
    /// iterations.
    /// </para>
    /// </summary>
    public static class ContactSolver2D
    {
        /// <summary>
        /// Prepares a manifold for the velocity iterations.
        /// </summary>
        /// <param name="manifold">The manifold to prepare.</param>
        /// <param name="restitutionThreshold">Approach speed below which a contact is inelastic.</param>
        /// <param name="warmStart">Apply the impulses remembered from the previous step.</param>
        public static void PreStep(Manifold2D manifold, double restitutionThreshold = 1.0, bool warmStart = true)
        {
            Prepare(manifold, restitutionThreshold);

            if (warmStart)
            {
                WarmStart(manifold);
                return;
            }

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                manifold.Contacts[i].NormalImpulse = 0;
                manifold.Contacts[i].TangentImpulse = 0;
            }
        }

        /// <summary>
        /// Computes the effective masses and restitution targets of a manifold without
        /// applying any impulse. The accumulated impulses are left untouched.
        /// </summary>
        public static void Prepare(Manifold2D manifold, double restitutionThreshold = 1.0)
        {
            var a = manifold.BodyA;
            var b = manifold.BodyB;

            double invMassA = a.SolverInverseMass, invInertiaA = a.SolverInverseInertia;
            double invMassB = b.SolverInverseMass, invInertiaB = b.SolverInverseInertia;
            double invMassSum = invMassA + invMassB;

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                ref var c = ref manifold.Contacts[i];

                c.RA = c.Point - a.Position;
                c.RB = c.Point - b.Position;

                var normal = c.Normal;
                var tangent = new Vector2D(-normal.Y, normal.X);

                double rnA = Vector2D.Cross(c.RA, normal);
                double rnB = Vector2D.Cross(c.RB, normal);
                double kNormal = invMassSum + invInertiaA * rnA * rnA + invInertiaB * rnB * rnB;
                c.NormalMass = kNormal > PhysicsConstants.Epsilon ? 1.0 / kNormal : 0;

                double rtA = Vector2D.Cross(c.RA, tangent);
                double rtB = Vector2D.Cross(c.RB, tangent);
                double kTangent = invMassSum + invInertiaA * rtA * rtA + invInertiaB * rtB * rtB;
                c.TangentMass = kTangent > PhysicsConstants.Epsilon ? 1.0 / kTangent : 0;

                // Restitution uses the approach speed before any impulse is applied.
                double approach = Vector2D.Dot(RelativeVelocity(a, b, c.RA, c.RB), normal);
                c.TargetNormalVelocity = approach < -restitutionThreshold ? -manifold.Restitution * approach : 0;
            }
        }

        /// <summary>
        /// Applies the impulses remembered from the previous step.
        /// <para>
        /// When several manifolds are solved together, prepare all of them first
        /// (<c>PreStep(..., warmStart: false)</c> clears the impulses, so use
        /// <see cref="Prepare"/>) and warm start afterwards: restitution is based on
        /// the approach speed, which must be measured before any impulse is applied.
        /// </para>
        /// </summary>
        public static void WarmStart(Manifold2D manifold)
        {
            var a = manifold.BodyA;
            var b = manifold.BodyB;

            double invMassA = a.SolverInverseMass, invInertiaA = a.SolverInverseInertia;
            double invMassB = b.SolverInverseMass, invInertiaB = b.SolverInverseInertia;

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                ref var c = ref manifold.Contacts[i];
                if (c.NormalImpulse == 0 && c.TangentImpulse == 0)
                    continue;

                var tangent = new Vector2D(-c.Normal.Y, c.Normal.X);
                var impulse = c.Normal * c.NormalImpulse + tangent * c.TangentImpulse;
                a.ApplySolverImpulse(-impulse, c.RA, invMassA, invInertiaA);
                b.ApplySolverImpulse(impulse, c.RB, invMassB, invInertiaB);
            }
        }

        /// <summary>
        /// One velocity iteration over the contacts of a manifold.
        /// </summary>
        public static void SolveVelocity(Manifold2D manifold)
        {
            var a = manifold.BodyA;
            var b = manifold.BodyB;

            double invMassA = a.SolverInverseMass, invInertiaA = a.SolverInverseInertia;
            double invMassB = b.SolverInverseMass, invInertiaB = b.SolverInverseInertia;
            if (invMassA + invMassB <= 0)
                return;

            double friction = manifold.Friction;

            // Friction first, so that non-penetration has the last word.
            for (int i = 0; i < manifold.ContactCount; i++)
            {
                ref var c = ref manifold.Contacts[i];
                var tangent = new Vector2D(-c.Normal.Y, c.Normal.X);

                double vt = Vector2D.Dot(RelativeVelocity(a, b, c.RA, c.RB), tangent);
                double lambda = -vt * c.TangentMass;

                // Coulomb: the accumulated friction impulse stays inside the cone.
                double limit = friction * c.NormalImpulse;
                double accumulated = Math.Clamp(c.TangentImpulse + lambda, -limit, limit);
                lambda = accumulated - c.TangentImpulse;
                c.TangentImpulse = accumulated;

                if (lambda != 0)
                {
                    var impulse = tangent * lambda;
                    a.ApplySolverImpulse(-impulse, c.RA, invMassA, invInertiaA);
                    b.ApplySolverImpulse(impulse, c.RB, invMassB, invInertiaB);
                }
            }

            if (manifold.ContactCount == 2 &&
                SolveNormalBlock(manifold, a, b, invMassA, invInertiaA, invMassB, invInertiaB))
                return;

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                ref var c = ref manifold.Contacts[i];

                double vn = Vector2D.Dot(RelativeVelocity(a, b, c.RA, c.RB), c.Normal);
                double lambda = -(vn - c.TargetNormalVelocity) * c.NormalMass;

                // Contacts can only push.
                double accumulated = Math.Max(c.NormalImpulse + lambda, 0);
                lambda = accumulated - c.NormalImpulse;
                c.NormalImpulse = accumulated;

                if (lambda != 0)
                {
                    var impulse = c.Normal * lambda;
                    a.ApplySolverImpulse(-impulse, c.RA, invMassA, invInertiaA);
                    b.ApplySolverImpulse(impulse, c.RB, invMassB, invInertiaB);
                }
            }
        }

        /// <summary>
        /// Solves both normal constraints of a two-point manifold together.
        /// <para>
        /// Solving the two points one after the other pushes on one corner first and
        /// leaves a small spin behind, which makes tall stacks sway. Treating them as a
        /// 2x2 complementarity problem (enumerating which points stay in contact) gives
        /// the exact answer for the pair. Returns false if the system is too
        /// ill-conditioned, in which case the caller falls back to the sequential solve.
        /// </para>
        /// </summary>
        private static bool SolveNormalBlock(Manifold2D manifold, RigidBody2D a, RigidBody2D b,
            double invMassA, double invInertiaA, double invMassB, double invInertiaB)
        {
            ref var c1 = ref manifold.Contacts[0];
            ref var c2 = ref manifold.Contacts[1];
            var normal = c1.Normal;

            double rn1A = Vector2D.Cross(c1.RA, normal), rn1B = Vector2D.Cross(c1.RB, normal);
            double rn2A = Vector2D.Cross(c2.RA, normal), rn2B = Vector2D.Cross(c2.RB, normal);
            double invMassSum = invMassA + invMassB;

            double k11 = invMassSum + invInertiaA * rn1A * rn1A + invInertiaB * rn1B * rn1B;
            double k22 = invMassSum + invInertiaA * rn2A * rn2A + invInertiaB * rn2B * rn2B;
            double k12 = invMassSum + invInertiaA * rn1A * rn2A + invInertiaB * rn1B * rn2B;
            double det = k11 * k22 - k12 * k12;

            // The two points must not be (nearly) redundant.
            if (k11 * k11 >= 1000.0 * det || det <= 0)
                return false;

            double a1 = c1.NormalImpulse, a2 = c2.NormalImpulse;

            double vn1 = Vector2D.Dot(RelativeVelocity(a, b, c1.RA, c1.RB), normal);
            double vn2 = Vector2D.Dot(RelativeVelocity(a, b, c2.RA, c2.RB), normal);

            // Velocity error with the effect of the accumulated impulses removed.
            double b1 = vn1 - c1.TargetNormalVelocity - (k11 * a1 + k12 * a2);
            double b2 = vn2 - c2.TargetNormalVelocity - (k12 * a1 + k22 * a2);

            double x1, x2;

            // Case 1: both points stay in contact.
            x1 = -(k22 * b1 - k12 * b2) / det;
            x2 = -(k11 * b2 - k12 * b1) / det;
            if (!(x1 >= 0 && x2 >= 0))
            {
                // Case 2: only point 1 pushes.
                x1 = -b1 / k11;
                x2 = 0;
                if (!(x1 >= 0 && k12 * x1 + b2 >= 0))
                {
                    // Case 3: only point 2 pushes.
                    x1 = 0;
                    x2 = -b2 / k22;
                    if (!(x2 >= 0 && k12 * x2 + b1 >= 0))
                    {
                        // Case 4: neither pushes.
                        if (!(b1 >= 0 && b2 >= 0))
                            return false;
                        x1 = 0;
                        x2 = 0;
                    }
                }
            }

            var p1 = normal * (x1 - a1);
            var p2 = normal * (x2 - a2);
            c1.NormalImpulse = x1;
            c2.NormalImpulse = x2;

            a.ApplySolverImpulse(-p1, c1.RA, invMassA, invInertiaA);
            a.ApplySolverImpulse(-p2, c2.RA, invMassA, invInertiaA);
            b.ApplySolverImpulse(p1, c1.RB, invMassB, invInertiaB);
            b.ApplySolverImpulse(p2, c2.RB, invMassB, invInertiaB);
            return true;
        }

        /// <summary>
        /// One positional iteration: accumulates the translation and rotation needed
        /// to remove the remaining overlap in the bodies' correction scratch state.
        /// </summary>
        internal static void SolvePosition(Manifold2D manifold, double slop, double factor, double maxCorrection)
        {
            var a = manifold.BodyA;
            var b = manifold.BodyB;

            double invMassA = a.SolverInverseMass, invInertiaA = a.SolverInverseInertia;
            double invMassB = b.SolverInverseMass, invInertiaB = b.SolverInverseInertia;
            if (invMassA + invMassB <= 0)
                return;

            for (int i = 0; i < manifold.ContactCount; i++)
            {
                ref var c = ref manifold.Contacts[i];

                var displacementA = a.CorrectionLinear + Vector2D.Cross(a.CorrectionAngular, c.RA);
                var displacementB = b.CorrectionLinear + Vector2D.Cross(b.CorrectionAngular, c.RB);
                double remaining = c.Depth - slop - Vector2D.Dot(displacementB - displacementA, c.Normal);
                if (remaining <= 0)
                    continue;

                double lambda = Math.Min(remaining * factor, maxCorrection) * c.NormalMass;
                var push = c.Normal * lambda;

                a.CorrectionLinear -= push * invMassA;
                a.CorrectionAngular -= Vector2D.Cross(c.RA, push) * invInertiaA;
                b.CorrectionLinear += push * invMassB;
                b.CorrectionAngular += Vector2D.Cross(c.RB, push) * invInertiaB;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vector2D RelativeVelocity(RigidBody2D a, RigidBody2D b, Vector2D rA, Vector2D rB)
        {
            var velA = a.Velocity + Vector2D.Cross(a.AngularVelocity, rA);
            var velB = b.Velocity + Vector2D.Cross(b.AngularVelocity, rB);
            return velB - velA;
        }
    }
}
