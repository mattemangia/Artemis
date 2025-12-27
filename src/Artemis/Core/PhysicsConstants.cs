using System;

namespace Artemis.Core
{
    /// <summary>
    /// Physical constants and configuration values used throughout the engine.
    /// </summary>
    public static class PhysicsConstants
    {
        #region Mathematical Constants

        /// <summary>
        /// Small value used for floating-point comparisons.
        /// </summary>
        public const double Epsilon = 1e-10;

        /// <summary>
        /// Pi constant.
        /// </summary>
        public const double Pi = Math.PI;

        /// <summary>
        /// Two times Pi (full circle in radians).
        /// </summary>
        public const double TwoPi = 2.0 * Math.PI;

        /// <summary>
        /// Degrees to radians conversion factor.
        /// </summary>
        public const double DegToRad = Math.PI / 180.0;

        /// <summary>
        /// Radians to degrees conversion factor.
        /// </summary>
        public const double RadToDeg = 180.0 / Math.PI;

        #endregion

        #region Physical Constants

        /// <summary>
        /// Standard Earth gravity in m/s².
        /// </summary>
        public const double EarthGravity = 9.80665;

        /// <summary>
        /// Moon gravity in m/s².
        /// </summary>
        public const double MoonGravity = 1.62;

        /// <summary>
        /// Mars gravity in m/s².
        /// </summary>
        public const double MarsGravity = 3.72076;

        /// <summary>
        /// Universal gravitational constant in N⋅m²/kg².
        /// </summary>
        public const double GravitationalConstant = 6.67430e-11;

        /// <summary>
        /// Alias for the universal gravitational constant.
        /// </summary>
        public const double G = GravitationalConstant;

        /// <summary>
        /// Speed of light in m/s.
        /// </summary>
        public const double SpeedOfLight = 299792458.0;

        /// <summary>
        /// Boltzmann constant in J/K.
        /// </summary>
        public const double BoltzmannConstant = 1.380649e-23;

        /// <summary>
        /// Avogadro's number.
        /// </summary>
        public const double AvogadroNumber = 6.02214076e23;

        #endregion

        #region Default Simulation Values

        /// <summary>
        /// Default fixed timestep for physics simulation (60 Hz).
        /// </summary>
        public const double DefaultFixedTimeStep = 1.0 / 60.0;

        /// <summary>
        /// Default maximum substeps per frame.
        /// </summary>
        public const int DefaultMaxSubSteps = 8;

        /// <summary>
        /// Default velocity sleep threshold.
        /// </summary>
        public const double DefaultSleepVelocityThreshold = 0.01;

        /// <summary>
        /// Default angular velocity sleep threshold.
        /// </summary>
        public const double DefaultSleepAngularVelocityThreshold = 0.01;

        /// <summary>
        /// Default time before a body can sleep.
        /// </summary>
        public const double DefaultSleepTimeThreshold = 0.5;

        /// <summary>
        /// Default linear damping coefficient.
        /// </summary>
        public const double DefaultLinearDamping = 0.01;

        /// <summary>
        /// Default angular damping coefficient.
        /// </summary>
        public const double DefaultAngularDamping = 0.05;

        #endregion

        #region Default Material Properties

        /// <summary>
        /// Default coefficient of restitution (bounciness).
        /// </summary>
        public const double DefaultRestitution = 0.3;

        /// <summary>
        /// Default static friction coefficient.
        /// </summary>
        public const double DefaultStaticFriction = 0.5;

        /// <summary>
        /// Default dynamic friction coefficient.
        /// </summary>
        public const double DefaultDynamicFriction = 0.3;

        /// <summary>
        /// Default elasticity (Young's modulus) in Pa.
        /// </summary>
        public const double DefaultElasticity = 2.1e11; // Steel

        /// <summary>
        /// Default Poisson's ratio.
        /// </summary>
        public const double DefaultPoissonRatio = 0.3;

        /// <summary>
        /// Default yield strength in Pa (for plasticity).
        /// </summary>
        public const double DefaultYieldStrength = 2.5e8; // Steel

        /// <summary>
        /// Default density in kg/m³.
        /// </summary>
        public const double DefaultDensity = 1000.0; // Water

        #endregion
    }
}
