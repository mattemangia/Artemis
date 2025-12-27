using System;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Simulates a vortex/whirlpool/tornado force that pulls objects in a spiral pattern.
    /// </summary>
    public class VortexForce : IForce
    {
        /// <inheritdoc/>
        public string Id { get; }

        /// <inheritdoc/>
        public bool Enabled { get; set; } = true;

        /// <summary>
        /// Gets or sets the center of the vortex.
        /// </summary>
        public Vector3D Center { get; set; }

        /// <summary>
        /// Gets or sets the vortex axis (direction of rotation axis).
        /// </summary>
        public Vector3D Axis { get; set; }

        /// <summary>
        /// Gets or sets the tangential (rotational) strength.
        /// </summary>
        public double TangentialStrength { get; set; }

        /// <summary>
        /// Gets or sets the radial (inward/outward) strength.
        /// Positive = inward (whirlpool), Negative = outward (explosion).
        /// </summary>
        public double RadialStrength { get; set; }

        /// <summary>
        /// Gets or sets the axial strength (up/down along axis).
        /// Positive = along axis direction, Negative = opposite.
        /// </summary>
        public double AxialStrength { get; set; }

        /// <summary>
        /// Gets or sets the inner radius where force is maximum.
        /// </summary>
        public double InnerRadius { get; set; }

        /// <summary>
        /// Gets or sets the outer radius where force reaches zero.
        /// </summary>
        public double OuterRadius { get; set; }

        /// <summary>
        /// Gets or sets the height of the vortex effect.
        /// </summary>
        public double Height { get; set; }

        /// <summary>
        /// Gets or sets the falloff exponent for force calculation.
        /// Higher values = sharper falloff.
        /// </summary>
        public double FalloffExponent { get; set; }

        /// <summary>
        /// Gets or sets whether rotation is clockwise when viewed from axis direction.
        /// </summary>
        public bool Clockwise { get; set; }

        /// <summary>
        /// Creates a new vortex force.
        /// </summary>
        public VortexForce(
            Vector3D center,
            Vector3D? axis = null,
            double tangentialStrength = 50.0,
            double radialStrength = 20.0,
            string? id = null)
        {
            Id = id ?? $"Vortex_{Guid.NewGuid():N}";
            Center = center;
            Axis = axis?.Normalized ?? Vector3D.Up;
            TangentialStrength = tangentialStrength;
            RadialStrength = radialStrength;
            AxialStrength = 0;
            InnerRadius = 0.5;
            OuterRadius = 10.0;
            Height = double.PositiveInfinity;
            FalloffExponent = 1.0;
            Clockwise = true;
        }

        /// <inheritdoc/>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled)
                return Vector3D.Zero;

            // Vector from center to position
            var toPosition = position - Center;

            // Project onto plane perpendicular to axis
            var axisNorm = Axis.Normalized;
            double heightAlongAxis = Vector3D.Dot(toPosition, axisNorm);

            // Check height bounds
            if (Math.Abs(heightAlongAxis) > Height / 2)
                return Vector3D.Zero;

            // Radial vector (in plane perpendicular to axis)
            var radial = toPosition - axisNorm * heightAlongAxis;
            double distance = radial.Magnitude;

            // Check radius bounds
            if (distance > OuterRadius || distance < PhysicsConstants.Epsilon)
                return Vector3D.Zero;

            // Calculate falloff
            double falloff;
            if (distance < InnerRadius)
            {
                falloff = 1.0;
            }
            else
            {
                double t = (distance - InnerRadius) / (OuterRadius - InnerRadius);
                falloff = Math.Pow(1.0 - t, FalloffExponent);
            }

            var radialDir = radial / distance;

            // Tangential direction (perpendicular to radial, in the plane)
            var tangentDir = Vector3D.Cross(axisNorm, radialDir);
            if (!Clockwise)
                tangentDir = -tangentDir;

            // Calculate force components
            var force = Vector3D.Zero;

            // Tangential force (rotation)
            force += tangentDir * TangentialStrength * falloff;

            // Radial force (inward/outward)
            force -= radialDir * RadialStrength * falloff;

            // Axial force (up/down)
            force += axisNorm * AxialStrength * falloff;

            return force;
        }

        #region Fluent API

        /// <summary>
        /// Sets the vortex center.
        /// </summary>
        public VortexForce AtPosition(Vector3D center)
        {
            Center = center;
            return this;
        }

        /// <summary>
        /// Sets the rotation axis.
        /// </summary>
        public VortexForce WithAxis(Vector3D axis)
        {
            Axis = axis.Normalized;
            return this;
        }

        /// <summary>
        /// Sets the tangential (rotational) strength.
        /// </summary>
        public VortexForce WithRotation(double strength)
        {
            TangentialStrength = strength;
            return this;
        }

        /// <summary>
        /// Sets the radial (pull) strength.
        /// </summary>
        public VortexForce WithPull(double strength)
        {
            RadialStrength = strength;
            return this;
        }

        /// <summary>
        /// Sets the axial (lift) strength.
        /// </summary>
        public VortexForce WithLift(double strength)
        {
            AxialStrength = strength;
            return this;
        }

        /// <summary>
        /// Sets the effective radius range.
        /// </summary>
        public VortexForce WithRadius(double inner, double outer)
        {
            InnerRadius = inner;
            OuterRadius = outer;
            return this;
        }

        /// <summary>
        /// Sets the vortex height.
        /// </summary>
        public VortexForce WithHeight(double height)
        {
            Height = height;
            return this;
        }

        /// <summary>
        /// Sets the rotation direction.
        /// </summary>
        public VortexForce RotatingClockwise(bool clockwise = true)
        {
            Clockwise = clockwise;
            return this;
        }

        #endregion

        #region Factory Methods

        /// <summary>
        /// Creates a whirlpool (water vortex pulling down).
        /// </summary>
        public static VortexForce Whirlpool(Vector3D center, double radius = 5.0)
            => new VortexForce(center, Vector3D.Up, 30, 40)
                .WithRadius(0.3, radius)
                .WithLift(-20);

        /// <summary>
        /// Creates a tornado.
        /// </summary>
        public static VortexForce Tornado(Vector3D center, double radius = 10.0)
            => new VortexForce(center, Vector3D.Up, 100, 50)
                .WithRadius(1.0, radius)
                .WithLift(80)
                .WithHeight(50);

        /// <summary>
        /// Creates a dust devil (small tornado).
        /// </summary>
        public static VortexForce DustDevil(Vector3D center)
            => new VortexForce(center, Vector3D.Up, 30, 15)
                .WithRadius(0.2, 3.0)
                .WithLift(20)
                .WithHeight(10);

        /// <summary>
        /// Creates a drain (sink).
        /// </summary>
        public static VortexForce Drain(Vector3D center, double radius = 2.0)
            => new VortexForce(center, Vector3D.Up, 20, 60)
                .WithRadius(0.1, radius)
                .WithLift(-50);

        /// <summary>
        /// Creates a fan/propeller effect.
        /// </summary>
        public static VortexForce Propeller(Vector3D center, Vector3D direction)
            => new VortexForce(center, direction, 40, 5)
                .WithRadius(0.5, 3.0)
                .WithLift(30);

        #endregion
    }
}
