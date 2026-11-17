using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Collision
{
    /// <summary>
    /// Sequential-impulse contact solver.
    /// <para>
    /// All contacts of a step are solved together over several iterations with
    /// accumulated, clamped impulses. This is what lets stacks and resting contacts
    /// settle: a single pass over the collisions (as done by
    /// <see cref="CollisionResolver"/>) only ever sees one pair at a time.
    /// </para>
    /// <para>
    /// Usage per step: integrate velocities, detect collisions, call
    /// <see cref="SolveVelocities"/>, integrate positions, then call
    /// <see cref="SolvePositions"/> to remove the remaining overlap.
    /// </para>
    /// </summary>
    public sealed class ContactSolver
    {
        #region Nested Types

        private struct BodyState
        {
            public IPhysicsBody Body;
            public Vector3D Position;
            public Vector3D Velocity;
            public Vector3D AngularVelocity;
            public double InverseMass;
            public Matrix3x3 InverseInertia;
            public Vector3D Correction;
            public Vector3D AngularCorrection;
            public bool Dirty;
        }

        private struct Constraint
        {
            public int A;
            public int B;
            public Vector3D Normal;
            public Vector3D Tangent1;
            public Vector3D Tangent2;
            public Vector3D RA;
            public Vector3D RB;
            public double NormalMass;
            public double TangentMass1;
            public double TangentMass2;
            public double TargetNormalVelocity;
            public double StaticFriction;
            public double DynamicFriction;
            public double Penetration;
            public double NormalImpulse;
            public double TangentImpulse1;
            public double TangentImpulse2;
            public Vector3D LocalA;
        }

        #endregion

        #region Fields

        private struct CachedContact
        {
            public Vector3D LocalA;
            public double NormalImpulse;
            public Vector3D FrictionImpulse;
        }

        private struct CacheRange
        {
            public int Start;
            public int Count;
        }

        private struct ManifoldRange
        {
            public IPhysicsBody A;
            public IPhysicsBody B;
            public int Start;
            public int Count;
        }

        private readonly Dictionary<IPhysicsBody, int> _bodyIndex = new();

        // Impulses of the previous step, keyed by body pair, used for warm starting.
        private Dictionary<(IPhysicsBody, IPhysicsBody), CacheRange> _cache = new();
        private Dictionary<(IPhysicsBody, IPhysicsBody), CacheRange> _nextCache = new();
        private CachedContact[] _cachePoints = new CachedContact[128];
        private CachedContact[] _nextCachePoints = new CachedContact[128];
        private readonly List<ManifoldRange> _manifolds = new();
        private BodyState[] _bodies = new BodyState[64];
        private int _bodyCount;
        private Constraint[] _constraints = new Constraint[128];
        private int _constraintCount;

        #endregion

        #region Properties

        /// <summary>
        /// Number of velocity iterations. More iterations give stiffer stacks.
        /// </summary>
        public int VelocityIterations { get; set; } = 12;

        /// <summary>
        /// Number of position (penetration recovery) iterations.
        /// </summary>
        public int PositionIterations { get; set; } = 4;

        /// <summary>
        /// Whether the impulses of the previous step are reused as the starting point
        /// of the next one. Persistent contacts (stacks, piles, heavy-on-light) converge
        /// far better with it. Requires the same solver instance to be used every step.
        /// </summary>
        public bool WarmStarting { get; set; } = true;

        /// <summary>
        /// Distance (m) within which a contact point is considered the same point as
        /// one from the previous step.
        /// </summary>
        public double ContactMatchDistance { get; set; } = 0.05;

        /// <summary>
        /// Approach speed below which a contact is treated as perfectly inelastic.
        /// Prevents resting contacts from jittering.
        /// </summary>
        public double RestitutionThreshold { get; set; } = 1.0;

        /// <summary>
        /// Overlap that is tolerated without correction, in metres.
        /// </summary>
        public double PenetrationSlop { get; set; } = 0.005;

        /// <summary>
        /// Fraction of the remaining overlap removed per position iteration.
        /// </summary>
        public double PositionCorrectionFactor { get; set; } = 0.3;

        /// <summary>
        /// Largest positional correction applied to a contact in one step, in metres.
        /// </summary>
        public double MaxPositionCorrection { get; set; } = 0.25;

        #endregion

        #region Velocity Solve

        /// <summary>
        /// Solves the velocity constraints for a set of collisions and writes the
        /// resulting velocities back to the bodies. Sleeping bodies and non-dynamic
        /// bodies are treated as immovable.
        /// </summary>
        public void SolveVelocities(IReadOnlyList<CollisionInfo> collisions)
        {
            _bodyIndex.Clear();
            _bodyCount = 0;
            _constraintCount = 0;
            _manifolds.Clear();

            for (int i = 0; i < collisions.Count; i++)
                AddCollision(collisions[i]);

            if (WarmStarting)
            {
                for (int i = 0; i < _constraintCount; i++)
                {
                    ref var c = ref _constraints[i];
                    if (c.NormalImpulse == 0 && c.TangentImpulse1 == 0 && c.TangentImpulse2 == 0)
                        continue;

                    var impulse = c.Normal * c.NormalImpulse + c.Tangent1 * c.TangentImpulse1 + c.Tangent2 * c.TangentImpulse2;
                    ApplyImpulse(ref _bodies[c.A], ref _bodies[c.B], c.RA, c.RB, impulse);
                }
            }

            // Each iteration solves friction for every contact and then the normal
            // constraints, so non-penetration has the last word. Keeping the two passes
            // separate means friction reacts to the velocities left after *all* contact
            // points of a manifold have been balanced, instead of to the transient spin
            // produced by pushing on one corner at a time.
            for (int iteration = 0; iteration < VelocityIterations; iteration++)
            {
                for (int i = 0; i < _constraintCount; i++)
                    SolveFriction(ref _constraints[i]);

                for (int i = 0; i < _constraintCount; i++)
                    SolveNormal(ref _constraints[i]);
            }

            StoreImpulses();

            for (int i = 0; i < _bodyCount; i++)
            {
                ref var state = ref _bodies[i];
                if (!state.Dirty)
                    continue;

                if (state.Body is RigidBody rb)
                {
                    rb.SetVelocities(state.Velocity, state.AngularVelocity);
                }
                else
                {
                    state.Body.Velocity = state.Velocity;
                    state.Body.AngularVelocity = state.AngularVelocity;
                }
            }
        }

        private void AddCollision(in CollisionInfo collision)
        {
            if (!collision.HasCollision)
                return;

            var bodyA = collision.BodyA;
            var bodyB = collision.BodyB;

            bool movableA = IsMovable(bodyA);
            bool movableB = IsMovable(bodyB);
            if (!movableA && !movableB)
                return;

            int a = GetBodyIndex(bodyA, movableA);
            int b = GetBodyIndex(bodyB, movableB);

            ref var stateA = ref _bodies[a];
            ref var stateB = ref _bodies[b];

            var normal = collision.Normal;
            BuildTangents(normal, out var tangent1, out var tangent2);

            var matA = bodyA.Material;
            var matB = bodyB.Material;
            double restitution = 0.5 * ((matA?.Restitution ?? 0) + (matB?.Restitution ?? 0));
            double staticFriction = Math.Sqrt(Math.Max(0, (matA?.StaticFriction ?? 0) * (matB?.StaticFriction ?? 0)));
            double dynamicFriction = Math.Sqrt(Math.Max(0, (matA?.DynamicFriction ?? 0) * (matB?.DynamicFriction ?? 0)));
            // Kinetic friction can never exceed static friction.
            dynamicFriction = Math.Min(dynamicFriction, staticFriction);

            int count = collision.ContactCount;

            bool hasCache = false;
            CacheRange cached = default;
            Transform transformA = default;
            if (WarmStarting)
            {
                transformA = bodyA.Transform;
                hasCache = _cache.TryGetValue((bodyA, bodyB), out cached);
                _manifolds.Add(new ManifoldRange { A = bodyA, B = bodyB, Start = _constraintCount, Count = count });
            }

            double matchDistanceSq = ContactMatchDistance * ContactMatchDistance;

            for (int i = 0; i < count; i++)
            {
                var point = collision.GetContact(i);

                var c = new Constraint
                {
                    A = a,
                    B = b,
                    Normal = normal,
                    Tangent1 = tangent1,
                    Tangent2 = tangent2,
                    RA = point.Position - stateA.Position,
                    RB = point.Position - stateB.Position,
                    StaticFriction = staticFriction,
                    DynamicFriction = dynamicFriction,
                    Penetration = point.Penetration
                };

                c.NormalMass = 1.0 / EffectiveMassInverse(stateA, stateB, c.RA, c.RB, normal);
                c.TangentMass1 = 1.0 / EffectiveMassInverse(stateA, stateB, c.RA, c.RB, tangent1);
                c.TangentMass2 = 1.0 / EffectiveMassInverse(stateA, stateB, c.RA, c.RB, tangent2);

                // Restitution is based on the approach speed before any impulse is applied.
                double approach = Vector3D.Dot(RelativeVelocity(stateA, stateB, c.RA, c.RB), normal);
                c.TargetNormalVelocity = approach < -RestitutionThreshold ? -restitution * approach : 0;

                if (WarmStarting)
                {
                    c.LocalA = transformA.InverseTransformPoint(point.Position);

                    if (hasCache)
                    {
                        int best = -1;
                        double bestDistance = matchDistanceSq;
                        for (int k = cached.Start; k < cached.Start + cached.Count; k++)
                        {
                            double d = (_cachePoints[k].LocalA - c.LocalA).MagnitudeSquared;
                            if (d < bestDistance)
                            {
                                bestDistance = d;
                                best = k;
                            }
                        }

                        if (best >= 0)
                        {
                            // Slightly under-relaxed to stay robust when the manifold shifts.
                            const double carryOver = 0.9;
                            c.NormalImpulse = _cachePoints[best].NormalImpulse * carryOver;
                            c.TangentImpulse1 = Vector3D.Dot(_cachePoints[best].FrictionImpulse, tangent1) * carryOver;
                            c.TangentImpulse2 = Vector3D.Dot(_cachePoints[best].FrictionImpulse, tangent2) * carryOver;
                        }
                    }
                }

                if (_constraintCount == _constraints.Length)
                    Array.Resize(ref _constraints, _constraints.Length * 2);
                _constraints[_constraintCount++] = c;
            }
        }

        private void SolveFriction(ref Constraint c)
        {
            ref var a = ref _bodies[c.A];
            ref var b = ref _bodies[c.B];

            // Coulomb cone on the accumulated tangential impulse.
            var relVel = RelativeVelocity(a, b, c.RA, c.RB);
            double old1 = c.TangentImpulse1;
            double old2 = c.TangentImpulse2;
            double t1 = old1 - Vector3D.Dot(relVel, c.Tangent1) * c.TangentMass1;
            double t2 = old2 - Vector3D.Dot(relVel, c.Tangent2) * c.TangentMass2;

            double magnitude = Math.Sqrt(t1 * t1 + t2 * t2);
            if (magnitude > c.StaticFriction * c.NormalImpulse)
            {
                // Sliding: kinetic friction opposes the motion.
                double limit = c.DynamicFriction * c.NormalImpulse;
                double scale = magnitude > PhysicsConstants.Epsilon ? limit / magnitude : 0;
                t1 *= scale;
                t2 *= scale;
            }

            c.TangentImpulse1 = t1;
            c.TangentImpulse2 = t2;

            var friction = c.Tangent1 * (t1 - old1) + c.Tangent2 * (t2 - old2);
            if (friction.MagnitudeSquared > 0)
                ApplyImpulse(ref a, ref b, c.RA, c.RB, friction);
        }

        private void SolveNormal(ref Constraint c)
        {
            ref var a = ref _bodies[c.A];
            ref var b = ref _bodies[c.B];

            // Non-penetration impulse: accumulated and clamped to push only.
            var relVel = RelativeVelocity(a, b, c.RA, c.RB);
            double vn = Vector3D.Dot(relVel, c.Normal);
            double lambda = -(vn - c.TargetNormalVelocity) * c.NormalMass;
            double newImpulse = Math.Max(c.NormalImpulse + lambda, 0);
            lambda = newImpulse - c.NormalImpulse;
            c.NormalImpulse = newImpulse;
            if (lambda != 0)
                ApplyImpulse(ref a, ref b, c.RA, c.RB, c.Normal * lambda);
        }

        private void StoreImpulses()
        {
            // Swap in a fresh cache; pairs that are no longer in contact drop out.
            (_cache, _nextCache) = (_nextCache, _cache);
            (_cachePoints, _nextCachePoints) = (_nextCachePoints, _cachePoints);
            _cache.Clear();

            if (!WarmStarting)
                return;

            if (_cachePoints.Length < _constraintCount)
                _cachePoints = new CachedContact[Math.Max(_constraintCount, _cachePoints.Length * 2)];

            foreach (var manifold in _manifolds)
            {
                for (int i = manifold.Start; i < manifold.Start + manifold.Count; i++)
                {
                    ref var c = ref _constraints[i];
                    _cachePoints[i] = new CachedContact
                    {
                        LocalA = c.LocalA,
                        NormalImpulse = c.NormalImpulse,
                        FrictionImpulse = c.Tangent1 * c.TangentImpulse1 + c.Tangent2 * c.TangentImpulse2
                    };
                }

                _cache[(manifold.A, manifold.B)] = new CacheRange { Start = manifold.Start, Count = manifold.Count };
            }
        }

        #endregion

        #region Position Solve

        /// <summary>
        /// Pushes overlapping bodies apart. Call after positions have been integrated
        /// with the velocities produced by <see cref="SolveVelocities"/>. The correction
        /// is purely positional (translation and rotation), so it does not add energy
        /// to the system the way a velocity bias would.
        /// </summary>
        public void SolvePositions()
        {
            if (_constraintCount == 0)
                return;

            for (int i = 0; i < _bodyCount; i++)
            {
                _bodies[i].Correction = Vector3D.Zero;
                _bodies[i].AngularCorrection = Vector3D.Zero;
            }

            for (int iteration = 0; iteration < PositionIterations; iteration++)
            {
                for (int i = 0; i < _constraintCount; i++)
                {
                    ref var c = ref _constraints[i];
                    ref var a = ref _bodies[c.A];
                    ref var b = ref _bodies[c.B];

                    if (a.InverseMass + b.InverseMass <= 0)
                        continue;

                    // Overlap left at this contact point after the corrections so far.
                    var displacementA = a.Correction + Vector3D.Cross(a.AngularCorrection, c.RA);
                    var displacementB = b.Correction + Vector3D.Cross(b.AngularCorrection, c.RB);
                    double remaining = c.Penetration - PenetrationSlop
                                       - Vector3D.Dot(displacementB - displacementA, c.Normal);
                    if (remaining <= 0)
                        continue;

                    double lambda = Math.Min(remaining * PositionCorrectionFactor, MaxPositionCorrection) * c.NormalMass;
                    var push = c.Normal * lambda;

                    if (a.InverseMass > 0)
                    {
                        a.Correction -= push * a.InverseMass;
                        a.AngularCorrection -= a.InverseInertia * Vector3D.Cross(c.RA, push);
                    }

                    if (b.InverseMass > 0)
                    {
                        b.Correction += push * b.InverseMass;
                        b.AngularCorrection += b.InverseInertia * Vector3D.Cross(c.RB, push);
                    }
                }
            }

            for (int i = 0; i < _bodyCount; i++)
            {
                ref var state = ref _bodies[i];
                if (state.Correction.MagnitudeSquared <= 0 && state.AngularCorrection.MagnitudeSquared <= 0)
                    continue;

                if (state.Body is RigidBody rb)
                    rb.ApplyPositionCorrection(state.Correction, state.AngularCorrection);
                else
                    state.Body.Position += state.Correction;
            }
        }

        #endregion

        #region Helpers

        /// <summary>
        /// Forgets the impulses remembered for warm starting (for example after
        /// teleporting bodies or resetting a simulation).
        /// </summary>
        public void ClearCache()
        {
            _cache.Clear();
            _nextCache.Clear();
        }

        private static bool IsMovable(IPhysicsBody body)
            => body.BodyType == BodyType.Dynamic && body.IsActive && !body.IsSleeping && body.InverseMass > 0;

        private int GetBodyIndex(IPhysicsBody body, bool movable)
        {
            if (_bodyIndex.TryGetValue(body, out int index))
                return index;

            if (_bodyCount == _bodies.Length)
                Array.Resize(ref _bodies, _bodies.Length * 2);

            index = _bodyCount++;
            _bodies[index] = new BodyState
            {
                Body = body,
                Position = body.Position,
                // Kinematic bodies keep their velocity so that they can push and carry
                // dynamic bodies; they just never receive impulses.
                Velocity = body.BodyType == BodyType.Static ? Vector3D.Zero : body.Velocity,
                AngularVelocity = body.BodyType == BodyType.Static ? Vector3D.Zero : body.AngularVelocity,
                InverseMass = movable ? body.InverseMass : 0,
                InverseInertia = movable && body is RigidBody rb ? rb.WorldInverseInertiaTensor : Matrix3x3.Zero
            };
            _bodyIndex[body] = index;
            return index;
        }

        private static Vector3D RelativeVelocity(in BodyState a, in BodyState b, Vector3D rA, Vector3D rB)
        {
            var velA = a.Velocity + Vector3D.Cross(a.AngularVelocity, rA);
            var velB = b.Velocity + Vector3D.Cross(b.AngularVelocity, rB);
            return velB - velA;
        }

        private static double EffectiveMassInverse(in BodyState a, in BodyState b, Vector3D rA, Vector3D rB, Vector3D direction)
        {
            double k = a.InverseMass + b.InverseMass;

            var rnA = Vector3D.Cross(rA, direction);
            k += Vector3D.Dot(rnA, a.InverseInertia * rnA);

            var rnB = Vector3D.Cross(rB, direction);
            k += Vector3D.Dot(rnB, b.InverseInertia * rnB);

            return Math.Max(k, PhysicsConstants.Epsilon);
        }

        private static void ApplyImpulse(ref BodyState a, ref BodyState b, Vector3D rA, Vector3D rB, Vector3D impulse)
        {
            if (a.InverseMass > 0)
            {
                a.Velocity -= impulse * a.InverseMass;
                a.AngularVelocity -= a.InverseInertia * Vector3D.Cross(rA, impulse);
                a.Dirty = true;
            }

            if (b.InverseMass > 0)
            {
                b.Velocity += impulse * b.InverseMass;
                b.AngularVelocity += b.InverseInertia * Vector3D.Cross(rB, impulse);
                b.Dirty = true;
            }
        }

        private static void BuildTangents(Vector3D normal, out Vector3D tangent1, out Vector3D tangent2)
        {
            // Pick the world axis least aligned with the normal to avoid a degenerate cross product.
            var reference = Math.Abs(normal.X) < 0.57735 ? Vector3D.Right : Vector3D.Up;
            tangent1 = Vector3D.Cross(normal, reference).Normalized;
            tangent2 = Vector3D.Cross(normal, tangent1);
        }

        #endregion
    }
}
