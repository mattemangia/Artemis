using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Collision
{
    /// <summary>
    /// Continuous collision detection helpers that keep fast bodies from tunnelling
    /// through thin obstacles between two discrete steps.
    /// </summary>
    public static class ContinuousCollision
    {
        /// <summary>
        /// Conservative advancement of the body's core sphere against the other
        /// bodies. Returns the part of <paramref name="dt"/> the body can travel before
        /// it would pass through something; the regular contact solver then handles
        /// the impact on the following step.
        /// </summary>
        public static double TimeOfImpact(RigidBody body, IReadOnlyList<IPhysicsBody> bodies, double dt)
        {
            if (!NeedsSweep(body, dt, out var swept))
                return dt;

            return TimeOfImpact(body, bodies, dt, swept);
        }

        /// <summary>
        /// Whether a body moves far enough in one step to skip past an obstacle, and
        /// if so the bounds of the volume its core sphere sweeps. Use the bounds to
        /// collect candidate bodies from a spatial structure and pass them to
        /// <see cref="TimeOfImpact(RigidBody, IReadOnlyList{IPhysicsBody}, double, AABB)"/>.
        /// </summary>
        public static bool NeedsSweep(RigidBody body, double dt, out AABB swept)
        {
            swept = default;

            double coreRadius = CoreRadius(body);
            var velocity = body.Velocity;

            // Slow bodies cannot skip past anything: discrete detection is enough.
            if (coreRadius <= 0 || velocity.MagnitudeSquared * dt * dt <= coreRadius * coreRadius)
                return false;

            var start = body.Position;
            var reach = new Vector3D(coreRadius);
            var end = start + velocity * dt;
            swept = new AABB(
                new Vector3D(Math.Min(start.X, end.X), Math.Min(start.Y, end.Y), Math.Min(start.Z, end.Z)) - reach,
                new Vector3D(Math.Max(start.X, end.X), Math.Max(start.Y, end.Y), Math.Max(start.Z, end.Z)) + reach);
            return true;
        }

        private static double CoreRadius(RigidBody body)
        {
            switch (body.ShapeType)
            {
                case CollisionShapeType.Sphere:
                    return body.Radius;
                case CollisionShapeType.Box:
                    var h = body.HalfExtents;
                    return Math.Min(h.X, Math.Min(h.Y, h.Z));
                default:
                    return 0;
            }
        }

        /// <summary>
        /// Time of impact of a fast body against a set of candidate bodies.
        /// </summary>
        /// <param name="body">The moving body.</param>
        /// <param name="bodies">Candidates; may contain duplicates and the body itself.</param>
        /// <param name="dt">Step duration.</param>
        /// <param name="swept">Swept bounds from <see cref="NeedsSweep"/>.</param>
        public static double TimeOfImpact(RigidBody body, IReadOnlyList<IPhysicsBody> bodies, double dt, AABB swept)
        {
            double coreRadius = CoreRadius(body);
            var velocity = body.Velocity;
            var start = body.Position;

            const double tolerance = 1e-4;
            double earliest = dt;

            for (int index = 0; index < bodies.Count; index++)
            {
                var other = bodies[index];
                if (ReferenceEquals(other, body) || !other.IsActive || other is not RigidBody otherRb)
                    continue;
                if (otherRb.ShapeType != CollisionShapeType.Sphere && otherRb.ShapeType != CollisionShapeType.Box)
                    continue;

                var otherVelocity = other.BodyType == BodyType.Static || other.IsSleeping
                    ? Vector3D.Zero
                    : other.Velocity;

                // Tunnelling needs a large *relative* motion. Bodies travelling together
                // (debris, a falling pile) are left to the discrete solver, which also
                // makes this the cheapest possible reject for the common case.
                var relativeVelocity = velocity - otherVelocity;
                double relativeSpeedSq = relativeVelocity.MagnitudeSquared;
                if (relativeSpeedSq * dt * dt <= coreRadius * coreRadius)
                    continue;

                var otherBounds = other.BoundingBox;
                if (otherVelocity.MagnitudeSquared > 0)
                {
                    var moved = otherVelocity * dt;
                    otherBounds.Encapsulate(otherBounds.Min + moved);
                    otherBounds.Encapsulate(otherBounds.Max + moved);
                }

                if (!swept.Intersects(otherBounds))
                    continue;

                double relativeSpeed = Math.Sqrt(relativeSpeedSq);

                double t = 0;
                for (int iteration = 0; iteration < 32 && t < earliest; iteration++)
                {
                    double distance = CollisionDetector.DistanceToBody(start + relativeVelocity * t, other) - coreRadius;
                    if (distance <= tolerance)
                    {
                        // Already touching at the start of the step is a job for the
                        // discrete solver; otherwise this is the time of impact.
                        if (iteration > 0)
                        {
                            // Stop just inside the surface so the overlap is detected next step.
                            double overshoot = Math.Min(0.25 * coreRadius, 0.01) / relativeSpeed;
                            earliest = Math.Min(earliest, t + overshoot);
                        }
                        break;
                    }

                    t += distance / relativeSpeed;
                }
            }

            return earliest;
        }
    }
}
