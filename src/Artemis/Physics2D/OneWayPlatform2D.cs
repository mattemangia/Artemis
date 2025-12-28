using System;
using System.Collections.Generic;

namespace Artemis.Physics2D
{
    /// <summary>
    /// One-way platform collision - objects can pass through from one direction.
    /// Commonly used for platforms you can jump up through.
    /// </summary>
    public class OneWayPlatformBehavior2D
    {
        public Vector2D AllowedDirection { get; set; } // Direction you CAN pass through
        public double Threshold { get; set; } = 0.1; // Velocity threshold

        public OneWayPlatformBehavior2D(Vector2D allowedDirection)
        {
            AllowedDirection = allowedDirection.Normalized;
        }

        /// <summary>
        /// Check if collision should be processed based on approach direction.
        /// </summary>
        public bool ShouldCollide(RigidBody2D moving, RigidBody2D platform, Vector2D collisionNormal)
        {
            // If moving in the allowed direction, pass through
            if (moving.Velocity.MagnitudeSquared < 0.0001)
                return true;

            double approachDot = Vector2D.Dot(moving.Velocity.Normalized, AllowedDirection);

            if (approachDot > Threshold)
            {
                // Moving in allowed direction - no collision
                return false;
            }

            // Check if approaching from the correct side
            double normalDot = Vector2D.Dot(collisionNormal, AllowedDirection);

            if (normalDot > 0)
            {
                // Collision normal points in allowed direction - allow pass through
                return false;
            }

            // Otherwise, collide normally
            return true;
        }

        /// <summary>
        /// Create a standard "jump-through" platform (pass through from below).
        /// </summary>
        public static OneWayPlatformBehavior2D CreateJumpThrough()
        {
            return new OneWayPlatformBehavior2D(new Vector2D(0, 1)); // Allow upward passage
        }

        /// <summary>
        /// Create a fall-through platform (pass through from above).
        /// </summary>
        public static OneWayPlatformBehavior2D CreateFallThrough()
        {
            return new OneWayPlatformBehavior2D(new Vector2D(0, -1)); // Allow downward passage
        }

        /// <summary>
        /// Create a left-to-right pass-through platform.
        /// </summary>
        public static OneWayPlatformBehavior2D CreateLeftToRight()
        {
            return new OneWayPlatformBehavior2D(new Vector2D(1, 0)); // Allow rightward passage
        }

        /// <summary>
        /// Create a right-to-left pass-through platform.
        /// </summary>
        public static OneWayPlatformBehavior2D CreateRightToLeft()
        {
            return new OneWayPlatformBehavior2D(new Vector2D(-1, 0)); // Allow leftward passage
        }
    }

    /// <summary>
    /// Extended RigidBody2D properties for one-way platforms.
    /// </summary>
    public static class OneWayPlatformExtensions2D
    {
        private static Dictionary<RigidBody2D, OneWayPlatformBehavior2D> _platformBehaviors = new();

        public static void SetOneWayPlatform(this RigidBody2D body, OneWayPlatformBehavior2D behavior)
        {
            _platformBehaviors[body] = behavior;
        }

        public static OneWayPlatformBehavior2D? GetOneWayPlatform(this RigidBody2D body)
        {
            return _platformBehaviors.TryGetValue(body, out var behavior) ? behavior : null;
        }

        public static bool IsOneWayPlatform(this RigidBody2D body)
        {
            return _platformBehaviors.ContainsKey(body);
        }

        public static void RemoveOneWayPlatform(this RigidBody2D body)
        {
            _platformBehaviors.Remove(body);
        }

        /// <summary>
        /// Check if collision between two bodies should be ignored due to one-way platform.
        /// </summary>
        public static bool ShouldIgnoreCollision(RigidBody2D bodyA, RigidBody2D bodyB, Vector2D normal)
        {
            var platformA = bodyA.GetOneWayPlatform();
            var platformB = bodyB.GetOneWayPlatform();

            if (platformA != null && !platformA.ShouldCollide(bodyB, bodyA, normal))
                return true;

            if (platformB != null && !platformB.ShouldCollide(bodyA, bodyB, -normal))
                return true;

            return false;
        }

        /// <summary>
        /// Make a body into a jump-through platform.
        /// </summary>
        public static void MakeJumpThrough(this RigidBody2D body)
        {
            body.SetOneWayPlatform(OneWayPlatformBehavior2D.CreateJumpThrough());
        }

        /// <summary>
        /// Make a body into a fall-through platform.
        /// </summary>
        public static void MakeFallThrough(this RigidBody2D body)
        {
            body.SetOneWayPlatform(OneWayPlatformBehavior2D.CreateFallThrough());
        }
    }

    /// <summary>
    /// Sleep system for 2D physics bodies.
    /// Puts inactive bodies to sleep to improve performance.
    /// </summary>
    public static class SleepingSystem2D
    {
        private const double SleepVelocityThreshold = 0.1;
        private const double SleepAngularVelocityThreshold = 0.05;
        private const double TimeToSleep = 0.5;

        private static Dictionary<RigidBody2D, double> _sleepTimers = new();

        public static void UpdateSleepState(RigidBody2D body, double deltaTime)
        {
            if (body.BodyType != BodyType2D.Dynamic || !body.CanSleep)
                return;

            double velocitySquared = body.Velocity.MagnitudeSquared;
            double angularVelocitySquared = body.AngularVelocity * body.AngularVelocity;

            if (velocitySquared < SleepVelocityThreshold * SleepVelocityThreshold &&
                angularVelocitySquared < SleepAngularVelocityThreshold * SleepAngularVelocityThreshold)
            {
                if (!_sleepTimers.ContainsKey(body))
                    _sleepTimers[body] = 0;

                _sleepTimers[body] += deltaTime;

                if (_sleepTimers[body] > TimeToSleep)
                {
                    body.IsSleeping = true;
                    body.Velocity = Vector2D.Zero;
                    body.AngularVelocity = 0;
                }
            }
            else
            {
                _sleepTimers[body] = 0;
                body.IsSleeping = false;
            }
        }

        public static void WakeUp(RigidBody2D body)
        {
            if (body.BodyType != BodyType2D.Static)
            {
                body.IsSleeping = false;
                _sleepTimers[body] = 0;
            }
        }

        public static void WakeUpNearby(PhysicsWorld2D world, RigidBody2D body, double radius = 5.0)
        {
            Vector2D pos = body.Position;
            double radiusSquared = radius * radius;

            foreach (var other in world.Bodies)
            {
                if (other == body || other.BodyType == BodyType2D.Static)
                    continue;

                double distSquared = (other.Position - pos).MagnitudeSquared;
                if (distSquared < radiusSquared)
                {
                    WakeUp(other);
                }
            }
        }

        public static void ClearTimers()
        {
            _sleepTimers.Clear();
        }
    }
}
