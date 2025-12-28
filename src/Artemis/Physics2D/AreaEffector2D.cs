using System;
using System.Runtime.CompilerServices;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Base class for area-based force effectors in 2D.
    /// </summary>
    public abstract class AreaEffector2D
    {
        public string Id { get; set; } = Guid.NewGuid().ToString();
        public bool Enabled { get; set; } = true;
        public Vector2D Position { get; set; }
        public double Radius { get; set; }

        /// <summary>
        /// Checks if this effector affects the given body.
        /// </summary>
        public virtual bool AffectsBody(RigidBody2D body)
        {
            double dist = (body.Position - Position).Magnitude;
            return dist <= Radius;
        }

        /// <summary>
        /// Applies force to the body.
        /// </summary>
        public abstract void ApplyForce(RigidBody2D body, double dt);
    }

    /// <summary>
    /// Radial force effector (attraction/repulsion from center).
    /// </summary>
    public class RadialForceEffector2D : AreaEffector2D
    {
        public double Strength { get; set; }
        public RadialFalloff2D Falloff { get; set; } = RadialFalloff2D.Linear;
        public bool AffectsByCharge { get; set; } = false;

        public RadialForceEffector2D(Vector2D center, double radius, double strength,
            RadialFalloff2D falloff = RadialFalloff2D.Linear)
        {
            Position = center;
            Radius = radius;
            Strength = strength;
            Falloff = falloff;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override void ApplyForce(RigidBody2D body, double dt)
        {
            var delta = body.Position - Position;
            double distance = delta.Magnitude;

            if (distance < 0.001 || distance > Radius)
                return;

            double falloffFactor = Falloff switch
            {
                RadialFalloff2D.None => 1.0,
                RadialFalloff2D.Linear => 1.0 - distance / Radius,
                RadialFalloff2D.Quadratic => 1.0 - (distance * distance) / (Radius * Radius),
                RadialFalloff2D.InverseSquare => (Radius * Radius) / (distance * distance),
                _ => 1.0
            };

            var direction = delta / distance;
            double forceMag = Strength * falloffFactor;

            // If affects by charge, use body's charge property
            if (AffectsByCharge && body.UserData is ICharged charged)
            {
                forceMag *= charged.Charge;
            }

            body.ApplyForce(direction * forceMag * dt);
        }
    }

    /// <summary>
    /// Directional force effector (uniform force in a direction).
    /// </summary>
    public class DirectionalForceEffector2D : AreaEffector2D
    {
        public Vector2D Direction { get; set; }
        public double Strength { get; set; }

        public DirectionalForceEffector2D(Vector2D position, double radius, Vector2D direction, double strength)
        {
            Position = position;
            Radius = radius;
            Direction = direction.Normalized;
            Strength = strength;
        }

        public override void ApplyForce(RigidBody2D body, double dt)
        {
            body.ApplyForce(Direction * Strength * dt);
        }
    }

    /// <summary>
    /// Vortex force effector (tangential force around center).
    /// </summary>
    public class VortexEffector2D : AreaEffector2D
    {
        public double Strength { get; set; }
        public bool Clockwise { get; set; } = true;

        public VortexEffector2D(Vector2D center, double radius, double strength, bool clockwise = true)
        {
            Position = center;
            Radius = radius;
            Strength = strength;
            Clockwise = clockwise;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override void ApplyForce(RigidBody2D body, double dt)
        {
            var delta = body.Position - Position;
            double distance = delta.Magnitude;

            if (distance < 0.001 || distance > Radius)
                return;

            // Tangent direction
            Vector2D tangent;
            if (Clockwise)
                tangent = new Vector2D(-delta.Y, delta.X) / distance;
            else
                tangent = new Vector2D(delta.Y, -delta.X) / distance;

            double falloff = 1.0 - distance / Radius;
            body.ApplyForce(tangent * Strength * falloff * dt);
        }
    }

    /// <summary>
    /// Magnetic field effector (Lorentz force on moving charged particles).
    /// </summary>
    public class MagneticFieldEffector2D : AreaEffector2D
    {
        public double FieldStrength { get; set; }
        public bool IsUniform { get; set; } = true;

        public MagneticFieldEffector2D(Vector2D center, double radius, double fieldStrength)
        {
            Position = center;
            Radius = radius;
            FieldStrength = fieldStrength;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public override void ApplyForce(RigidBody2D body, double dt)
        {
            // Lorentz force: F = qv × B
            // In 2D, B is perpendicular to plane, so F is perpendicular to v
            var velocity = body.Velocity;
            if (velocity.MagnitudeSquared < 0.001)
                return;

            double charge = 1.0;
            if (body.UserData is ICharged charged)
                charge = charged.Charge;

            // Force perpendicular to velocity (circular motion)
            var force = new Vector2D(-velocity.Y, velocity.X) * FieldStrength * charge;

            // Apply centripetal adjustment for orbit around center
            if (!IsUniform)
            {
                var toCenter = Position - body.Position;
                double dist = toCenter.Magnitude;
                if (dist > 0.001 && dist < Radius)
                {
                    double centripetal = velocity.MagnitudeSquared / dist;
                    force += toCenter.Normalized * centripetal * body.Mass * 0.1;
                }
            }

            body.ApplyForce(force * dt);
        }
    }

    /// <summary>
    /// Buoyancy effector for fluid simulation.
    /// </summary>
    public class BuoyancyEffector2D : AreaEffector2D
    {
        public double FluidDensity { get; set; } = 1000.0;
        public double SurfaceY { get; set; }
        public double Drag { get; set; } = 0.5;

        public BuoyancyEffector2D(Vector2D position, double radius, double surfaceY, double fluidDensity = 1000.0)
        {
            Position = position;
            Radius = radius;
            SurfaceY = surfaceY;
            FluidDensity = fluidDensity;
        }

        public override bool AffectsBody(RigidBody2D body)
        {
            return body.Position.Y < SurfaceY && base.AffectsBody(body);
        }

        public override void ApplyForce(RigidBody2D body, double dt)
        {
            // Simplified buoyancy: F = ρVg
            double submergedDepth = Math.Max(0, SurfaceY - body.Position.Y);
            double buoyancy = FluidDensity * submergedDepth * 9.81 * 0.1; // Simplified

            body.ApplyForce(new Vector2D(0, buoyancy) * dt);

            // Drag
            body.Velocity *= (1.0 - Drag * dt);
        }
    }

    public enum RadialFalloff2D
    {
        None,
        Linear,
        Quadratic,
        InverseSquare
    }

    /// <summary>
    /// Interface for charged particles.
    /// </summary>
    public interface ICharged
    {
        double Charge { get; }
    }
}
