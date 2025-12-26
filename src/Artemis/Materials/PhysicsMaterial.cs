using System;
using Artemis.Core;

namespace Artemis.Materials
{
    /// <summary>
    /// Defines the physical properties of a material including elasticity, plasticity, and friction.
    /// </summary>
    public class PhysicsMaterial
    {
        #region Properties

        /// <summary>
        /// Gets or sets the name of the material.
        /// </summary>
        public string Name { get; set; }

        /// <summary>
        /// Gets or sets the density in kg/m³.
        /// </summary>
        public double Density { get; set; }

        /// <summary>
        /// Gets or sets the coefficient of restitution (bounciness).
        /// 0 = perfectly inelastic, 1 = perfectly elastic.
        /// </summary>
        public double Restitution { get; set; }

        /// <summary>
        /// Gets or sets the static friction coefficient.
        /// </summary>
        public double StaticFriction { get; set; }

        /// <summary>
        /// Gets or sets the dynamic (kinetic) friction coefficient.
        /// </summary>
        public double DynamicFriction { get; set; }

        /// <summary>
        /// Gets or sets Young's modulus (elasticity) in Pa.
        /// Higher values = stiffer material.
        /// </summary>
        public double YoungsModulus { get; set; }

        /// <summary>
        /// Gets or sets Poisson's ratio (lateral strain / axial strain).
        /// Typically between 0 and 0.5.
        /// </summary>
        public double PoissonRatio { get; set; }

        /// <summary>
        /// Gets or sets the yield strength in Pa (for plasticity).
        /// Stress above this value causes permanent deformation.
        /// </summary>
        public double YieldStrength { get; set; }

        /// <summary>
        /// Gets or sets the ultimate tensile strength in Pa.
        /// Stress above this value causes fracture.
        /// </summary>
        public double UltimateTensileStrength { get; set; }

        /// <summary>
        /// Gets or sets the hardening coefficient for plastic deformation.
        /// </summary>
        public double HardeningCoefficient { get; set; }

        /// <summary>
        /// Gets or sets the ductility (ability to deform under tensile stress).
        /// 0 = brittle, 1 = very ductile.
        /// </summary>
        public double Ductility { get; set; }

        /// <summary>
        /// Gets or sets the thermal conductivity in W/(m·K).
        /// </summary>
        public double ThermalConductivity { get; set; }

        /// <summary>
        /// Gets or sets the specific heat capacity in J/(kg·K).
        /// </summary>
        public double SpecificHeat { get; set; }

        /// <summary>
        /// Gets or sets the melting point in Kelvin.
        /// </summary>
        public double MeltingPoint { get; set; }

        #endregion

        #region Derived Properties

        /// <summary>
        /// Gets the shear modulus (rigidity) in Pa.
        /// G = E / (2 * (1 + ν))
        /// </summary>
        public double ShearModulus => YoungsModulus / (2.0 * (1.0 + PoissonRatio));

        /// <summary>
        /// Gets the bulk modulus in Pa.
        /// K = E / (3 * (1 - 2ν))
        /// </summary>
        public double BulkModulus
        {
            get
            {
                double denom = 3.0 * (1.0 - 2.0 * PoissonRatio);
                return Math.Abs(denom) > PhysicsConstants.Epsilon
                    ? YoungsModulus / denom
                    : double.PositiveInfinity;
            }
        }

        /// <summary>
        /// Gets the Lamé's first parameter.
        /// λ = νE / ((1 + ν)(1 - 2ν))
        /// </summary>
        public double LameFirst
        {
            get
            {
                double denom = (1.0 + PoissonRatio) * (1.0 - 2.0 * PoissonRatio);
                return Math.Abs(denom) > PhysicsConstants.Epsilon
                    ? PoissonRatio * YoungsModulus / denom
                    : 0;
            }
        }

        /// <summary>
        /// Gets the Lamé's second parameter (same as shear modulus).
        /// </summary>
        public double LameSecond => ShearModulus;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new physics material with default properties.
        /// </summary>
        public PhysicsMaterial(string name = "Default")
        {
            Name = name;
            Density = PhysicsConstants.DefaultDensity;
            Restitution = PhysicsConstants.DefaultRestitution;
            StaticFriction = PhysicsConstants.DefaultStaticFriction;
            DynamicFriction = PhysicsConstants.DefaultDynamicFriction;
            YoungsModulus = PhysicsConstants.DefaultElasticity;
            PoissonRatio = PhysicsConstants.DefaultPoissonRatio;
            YieldStrength = PhysicsConstants.DefaultYieldStrength;
            UltimateTensileStrength = PhysicsConstants.DefaultYieldStrength * 1.5;
            HardeningCoefficient = 0.1;
            Ductility = 0.5;
            ThermalConductivity = 50.0;
            SpecificHeat = 500.0;
            MeltingPoint = 1500.0;
        }

        /// <summary>
        /// Creates a new physics material with specified basic properties.
        /// </summary>
        public PhysicsMaterial(
            string name,
            double density,
            double restitution = 0.3,
            double staticFriction = 0.5,
            double dynamicFriction = 0.3)
        {
            Name = name;
            Density = density;
            Restitution = restitution;
            StaticFriction = staticFriction;
            DynamicFriction = dynamicFriction;
            YoungsModulus = PhysicsConstants.DefaultElasticity;
            PoissonRatio = PhysicsConstants.DefaultPoissonRatio;
            YieldStrength = PhysicsConstants.DefaultYieldStrength;
            UltimateTensileStrength = PhysicsConstants.DefaultYieldStrength * 1.5;
            HardeningCoefficient = 0.1;
            Ductility = 0.5;
            ThermalConductivity = 50.0;
            SpecificHeat = 500.0;
            MeltingPoint = 1500.0;
        }

        #endregion

        #region Methods

        /// <summary>
        /// Calculates the stress response for a given strain.
        /// Includes elastic and plastic behavior.
        /// </summary>
        /// <param name="strain">The strain (deformation ratio).</param>
        /// <param name="plasticStrain">Current accumulated plastic strain.</param>
        /// <returns>The stress in Pa and updated plastic strain.</returns>
        public (double stress, double newPlasticStrain) CalculateStress(double strain, double plasticStrain)
        {
            // Elastic strain
            double elasticStrain = strain - plasticStrain;
            double trialStress = YoungsModulus * elasticStrain;

            // Check for plastic yielding
            double yieldStress = YieldStrength + HardeningCoefficient * plasticStrain;

            if (Math.Abs(trialStress) <= yieldStress)
            {
                // Elastic response
                return (trialStress, plasticStrain);
            }

            // Plastic response
            double sign = Math.Sign(trialStress);
            double stress = sign * yieldStress;
            double newPlasticStrain = strain - stress / YoungsModulus;

            return (stress, newPlasticStrain);
        }

        /// <summary>
        /// Checks if the material would fracture at the given stress.
        /// </summary>
        public bool WouldFracture(double stress)
            => Math.Abs(stress) >= UltimateTensileStrength;

        /// <summary>
        /// Combines two materials for collision response.
        /// </summary>
        public static PhysicsMaterial Combine(PhysicsMaterial a, PhysicsMaterial b)
        {
            return new PhysicsMaterial("Combined")
            {
                Density = (a.Density + b.Density) * 0.5,
                Restitution = Math.Sqrt(a.Restitution * b.Restitution), // Geometric mean
                StaticFriction = Math.Sqrt(a.StaticFriction * b.StaticFriction),
                DynamicFriction = Math.Sqrt(a.DynamicFriction * b.DynamicFriction),
                YoungsModulus = 2.0 * a.YoungsModulus * b.YoungsModulus /
                               (a.YoungsModulus + b.YoungsModulus), // Harmonic mean
                PoissonRatio = (a.PoissonRatio + b.PoissonRatio) * 0.5,
                YieldStrength = Math.Min(a.YieldStrength, b.YieldStrength),
                UltimateTensileStrength = Math.Min(a.UltimateTensileStrength, b.UltimateTensileStrength),
                Ductility = (a.Ductility + b.Ductility) * 0.5
            };
        }

        /// <summary>
        /// Creates a copy of this material.
        /// </summary>
        public PhysicsMaterial Clone()
        {
            return new PhysicsMaterial(Name)
            {
                Density = Density,
                Restitution = Restitution,
                StaticFriction = StaticFriction,
                DynamicFriction = DynamicFriction,
                YoungsModulus = YoungsModulus,
                PoissonRatio = PoissonRatio,
                YieldStrength = YieldStrength,
                UltimateTensileStrength = UltimateTensileStrength,
                HardeningCoefficient = HardeningCoefficient,
                Ductility = Ductility,
                ThermalConductivity = ThermalConductivity,
                SpecificHeat = SpecificHeat,
                MeltingPoint = MeltingPoint
            };
        }

        public override string ToString() => $"Material({Name})";

        #endregion
    }
}
