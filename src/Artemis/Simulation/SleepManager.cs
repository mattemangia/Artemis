using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Collision;
using Artemis.Core;

namespace Artemis.Simulation
{
    /// <summary>
    /// Island based sleeping.
    /// <para>
    /// Dynamic bodies that touch each other form an island. An island is put to sleep
    /// only when all of its bodies have been at rest for a while, and it is woken as a
    /// whole as soon as something moving touches it. Deciding per body instead would
    /// freeze the bottom of a stack while the top is still settling, or leave a
    /// sleeping body hanging in the air after its support was knocked away.
    /// </para>
    /// </summary>
    public sealed class SleepManager
    {
        private readonly Dictionary<IPhysicsBody, int> _islandIndex = new();
        private int[] _islandParent = new int[64];
        private byte[] _islandFlag = new byte[64];

        /// <summary>
        /// Linear speed (m/s) below which a body counts as resting.
        /// </summary>
        public double SleepVelocityThreshold { get; set; } = 0.1;

        /// <summary>
        /// Angular speed (rad/s) below which a body counts as resting.
        /// </summary>
        public double SleepAngularVelocityThreshold { get; set; } = 0.25;

        /// <summary>
        /// Time (s) a whole island has to rest before it is put to sleep.
        /// </summary>
        public double SleepTimeThreshold { get; set; } = PhysicsConstants.DefaultSleepTimeThreshold;

        /// <summary>
        /// Whether a body is a kinematic body that is currently moving.
        /// </summary>
        public static bool IsMovingKinematic(IPhysicsBody body)
            => body.BodyType == BodyType.Kinematic &&
               (body.Velocity.MagnitudeSquared > PhysicsConstants.Epsilon ||
                body.AngularVelocity.MagnitudeSquared > PhysicsConstants.Epsilon);

        private int FindRoot(int index)
        {
            while (_islandParent[index] != index)
            {
                _islandParent[index] = _islandParent[_islandParent[index]];
                index = _islandParent[index];
            }
            return index;
        }

        /// <summary>
        /// Groups dynamic bodies that touch each other. Must be called once per step,
        /// after collision detection and before the other methods.
        /// </summary>
        public void BuildIslands(IReadOnlyList<IPhysicsBody> bodies, IReadOnlyList<CollisionInfo> collisions)
        {
            int count = bodies.Count;
            if (_islandParent.Length < count)
            {
                _islandParent = new int[Math.Max(count, _islandParent.Length * 2)];
                _islandFlag = new byte[_islandParent.Length];
            }

            _islandIndex.Clear();
            for (int i = 0; i < count; i++)
            {
                _islandParent[i] = i;
                _islandIndex[bodies[i]] = i;
            }

            foreach (var collision in collisions)
            {
                if (collision.BodyA.BodyType != BodyType.Dynamic || collision.BodyB.BodyType != BodyType.Dynamic)
                    continue;

                int a = FindRoot(_islandIndex[collision.BodyA]);
                int b = FindRoot(_islandIndex[collision.BodyB]);
                if (a != b)
                    _islandParent[b] = a;
            }
        }

        /// <summary>
        /// Wakes every sleeping body whose island contains a moving body or is touched
        /// by a moving kinematic body. Call before solving contacts.
        /// </summary>
        public void WakeDisturbedIslands(IReadOnlyList<IPhysicsBody> bodies, IReadOnlyList<CollisionInfo> collisions)
        {
            int count = bodies.Count;
            Array.Clear(_islandFlag, 0, count);
            bool any = false;

            // An island is disturbed when one of its awake members is actually moving...
            for (int i = 0; i < count; i++)
            {
                var body = bodies[i];
                if (body.BodyType != BodyType.Dynamic || !body.IsActive || body.IsSleeping)
                    continue;

                if (!IsResting(body))
                {
                    _islandFlag[FindRoot(i)] = 1;
                    any = true;
                }
            }

            // ...or when a moving kinematic body touches it.
            foreach (var collision in collisions)
            {
                IPhysicsBody? sleeper = null;
                if (IsMovingKinematic(collision.BodyA)) sleeper = collision.BodyB;
                else if (IsMovingKinematic(collision.BodyB)) sleeper = collision.BodyA;

                if (sleeper != null && sleeper.BodyType == BodyType.Dynamic)
                {
                    _islandFlag[FindRoot(_islandIndex[sleeper])] = 1;
                    any = true;
                }
            }

            if (!any)
                return;

            for (int i = 0; i < count; i++)
            {
                var body = bodies[i];
                if (body.IsSleeping && body.BodyType == BodyType.Dynamic && _islandFlag[FindRoot(i)] == 1)
                    body.WakeUp();
            }
        }

        private bool IsResting(IPhysicsBody body)
            => body.Velocity.MagnitudeSquared < SleepVelocityThreshold * SleepVelocityThreshold &&
               body.AngularVelocity.MagnitudeSquared < SleepAngularVelocityThreshold * SleepAngularVelocityThreshold;

        /// <summary>
        /// Advances the sleep timers and puts islands to sleep once every body in them
        /// has been at rest for <see cref="SleepTimeThreshold"/>. Call at the end of a step.
        /// </summary>
        public void UpdateSleeping(IReadOnlyList<IPhysicsBody> bodies, double dt)
        {
            int count = bodies.Count;
            Array.Clear(_islandFlag, 0, count);

            // Flag every island that contains a body which is not ready to sleep.
            for (int i = 0; i < count; i++)
            {
                if (bodies[i] is not RigidBody rb || rb.BodyType != BodyType.Dynamic || !rb.IsActive || rb.IsSleeping)
                    continue;

                if (rb.CanSleep && IsResting(rb))
                    rb.SleepTimer += dt;
                else
                    rb.SleepTimer = 0;

                if (rb.SleepTimer < SleepTimeThreshold)
                    _islandFlag[FindRoot(i)] = 1;
            }

            for (int i = 0; i < count; i++)
            {
                if (bodies[i] is RigidBody rb && rb.BodyType == BodyType.Dynamic && rb.IsActive && !rb.IsSleeping &&
                    _islandFlag[FindRoot(i)] == 0)
                {
                    rb.Sleep();
                }
            }
        }
    }
}
