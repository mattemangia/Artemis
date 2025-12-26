namespace Artemis.Materials
{
    /// <summary>
    /// Provides preset material configurations for common physical materials.
    /// </summary>
    public static class MaterialPresets
    {
        #region Metals

        /// <summary>
        /// Creates a steel material.
        /// </summary>
        public static PhysicsMaterial Steel() => new("Steel")
        {
            Density = 7850,
            Restitution = 0.6,
            StaticFriction = 0.74,
            DynamicFriction = 0.57,
            YoungsModulus = 2.0e11,
            PoissonRatio = 0.30,
            YieldStrength = 2.5e8,
            UltimateTensileStrength = 4.0e8,
            HardeningCoefficient = 0.15,
            Ductility = 0.7,
            ThermalConductivity = 50.2,
            SpecificHeat = 490,
            MeltingPoint = 1643
        };

        /// <summary>
        /// Creates an aluminum material.
        /// </summary>
        public static PhysicsMaterial Aluminum() => new("Aluminum")
        {
            Density = 2700,
            Restitution = 0.5,
            StaticFriction = 0.61,
            DynamicFriction = 0.47,
            YoungsModulus = 6.9e10,
            PoissonRatio = 0.33,
            YieldStrength = 2.75e8,
            UltimateTensileStrength = 3.1e8,
            HardeningCoefficient = 0.1,
            Ductility = 0.8,
            ThermalConductivity = 237,
            SpecificHeat = 897,
            MeltingPoint = 933
        };

        /// <summary>
        /// Creates a copper material.
        /// </summary>
        public static PhysicsMaterial Copper() => new("Copper")
        {
            Density = 8960,
            Restitution = 0.45,
            StaticFriction = 0.53,
            DynamicFriction = 0.36,
            YoungsModulus = 1.1e11,
            PoissonRatio = 0.34,
            YieldStrength = 7.0e7,
            UltimateTensileStrength = 2.2e8,
            HardeningCoefficient = 0.12,
            Ductility = 0.85,
            ThermalConductivity = 401,
            SpecificHeat = 385,
            MeltingPoint = 1358
        };

        /// <summary>
        /// Creates an iron material.
        /// </summary>
        public static PhysicsMaterial Iron() => new("Iron")
        {
            Density = 7874,
            Restitution = 0.55,
            StaticFriction = 0.7,
            DynamicFriction = 0.4,
            YoungsModulus = 2.11e11,
            PoissonRatio = 0.29,
            YieldStrength = 8.0e7,
            UltimateTensileStrength = 3.5e8,
            HardeningCoefficient = 0.1,
            Ductility = 0.6,
            ThermalConductivity = 80.4,
            SpecificHeat = 449,
            MeltingPoint = 1811
        };

        /// <summary>
        /// Creates a gold material.
        /// </summary>
        public static PhysicsMaterial Gold() => new("Gold")
        {
            Density = 19320,
            Restitution = 0.3,
            StaticFriction = 0.49,
            DynamicFriction = 0.4,
            YoungsModulus = 7.9e10,
            PoissonRatio = 0.44,
            YieldStrength = 2.0e7,
            UltimateTensileStrength = 1.2e8,
            HardeningCoefficient = 0.05,
            Ductility = 0.95,
            ThermalConductivity = 318,
            SpecificHeat = 129,
            MeltingPoint = 1337
        };

        #endregion

        #region Non-Metals

        /// <summary>
        /// Creates a rubber material.
        /// </summary>
        public static PhysicsMaterial Rubber() => new("Rubber")
        {
            Density = 1100,
            Restitution = 0.85,
            StaticFriction = 1.0,
            DynamicFriction = 0.8,
            YoungsModulus = 1.0e7,
            PoissonRatio = 0.49,
            YieldStrength = 1.5e7,
            UltimateTensileStrength = 2.5e7,
            HardeningCoefficient = 0.0,
            Ductility = 0.95,
            ThermalConductivity = 0.13,
            SpecificHeat = 2000,
            MeltingPoint = 453
        };

        /// <summary>
        /// Creates a wood material.
        /// </summary>
        public static PhysicsMaterial Wood() => new("Wood")
        {
            Density = 700,
            Restitution = 0.4,
            StaticFriction = 0.5,
            DynamicFriction = 0.3,
            YoungsModulus = 1.1e10,
            PoissonRatio = 0.35,
            YieldStrength = 4.0e7,
            UltimateTensileStrength = 8.0e7,
            HardeningCoefficient = 0.0,
            Ductility = 0.1,
            ThermalConductivity = 0.12,
            SpecificHeat = 1700,
            MeltingPoint = 573 // Ignition temperature
        };

        /// <summary>
        /// Creates a glass material.
        /// </summary>
        public static PhysicsMaterial Glass() => new("Glass")
        {
            Density = 2500,
            Restitution = 0.5,
            StaticFriction = 0.9,
            DynamicFriction = 0.4,
            YoungsModulus = 7.0e10,
            PoissonRatio = 0.22,
            YieldStrength = 3.3e7, // Glass doesn't really yield
            UltimateTensileStrength = 3.3e7,
            HardeningCoefficient = 0.0,
            Ductility = 0.0, // Brittle
            ThermalConductivity = 1.0,
            SpecificHeat = 840,
            MeltingPoint = 1673
        };

        /// <summary>
        /// Creates a concrete material.
        /// </summary>
        public static PhysicsMaterial Concrete() => new("Concrete")
        {
            Density = 2400,
            Restitution = 0.2,
            StaticFriction = 0.6,
            DynamicFriction = 0.45,
            YoungsModulus = 3.0e10,
            PoissonRatio = 0.2,
            YieldStrength = 3.0e7,
            UltimateTensileStrength = 4.0e7,
            HardeningCoefficient = 0.0,
            Ductility = 0.05,
            ThermalConductivity = 1.7,
            SpecificHeat = 880,
            MeltingPoint = 1773
        };

        /// <summary>
        /// Creates an ice material.
        /// </summary>
        public static PhysicsMaterial Ice() => new("Ice")
        {
            Density = 917,
            Restitution = 0.3,
            StaticFriction = 0.1,
            DynamicFriction = 0.03,
            YoungsModulus = 9.3e9,
            PoissonRatio = 0.33,
            YieldStrength = 1.0e6,
            UltimateTensileStrength = 2.0e6,
            HardeningCoefficient = 0.0,
            Ductility = 0.0,
            ThermalConductivity = 2.22,
            SpecificHeat = 2090,
            MeltingPoint = 273
        };

        #endregion

        #region Granular Materials

        /// <summary>
        /// Creates a sand material.
        /// </summary>
        public static PhysicsMaterial Sand() => new("Sand")
        {
            Density = 1600,
            Restitution = 0.1,
            StaticFriction = 0.6,
            DynamicFriction = 0.4,
            YoungsModulus = 1.0e8,
            PoissonRatio = 0.3,
            YieldStrength = 0, // Granular, no yield
            UltimateTensileStrength = 0,
            HardeningCoefficient = 0.0,
            Ductility = 0.0,
            ThermalConductivity = 0.25,
            SpecificHeat = 830,
            MeltingPoint = 1973 // Silica melting point
        };

        /// <summary>
        /// Creates a gravel material.
        /// </summary>
        public static PhysicsMaterial Gravel() => new("Gravel")
        {
            Density = 1800,
            Restitution = 0.15,
            StaticFriction = 0.7,
            DynamicFriction = 0.5,
            YoungsModulus = 5.0e8,
            PoissonRatio = 0.3,
            YieldStrength = 0,
            UltimateTensileStrength = 0,
            HardeningCoefficient = 0.0,
            Ductility = 0.0,
            ThermalConductivity = 0.5,
            SpecificHeat = 800,
            MeltingPoint = 1973
        };

        #endregion

        #region Soft Bodies

        /// <summary>
        /// Creates a soft tissue material (for biological simulations).
        /// </summary>
        public static PhysicsMaterial SoftTissue() => new("SoftTissue")
        {
            Density = 1050,
            Restitution = 0.3,
            StaticFriction = 0.4,
            DynamicFriction = 0.3,
            YoungsModulus = 1.0e5,
            PoissonRatio = 0.45,
            YieldStrength = 1.0e4,
            UltimateTensileStrength = 3.0e4,
            HardeningCoefficient = 0.0,
            Ductility = 0.9,
            ThermalConductivity = 0.5,
            SpecificHeat = 3500,
            MeltingPoint = 373
        };

        /// <summary>
        /// Creates a foam material.
        /// </summary>
        public static PhysicsMaterial Foam() => new("Foam")
        {
            Density = 100,
            Restitution = 0.7,
            StaticFriction = 0.5,
            DynamicFriction = 0.4,
            YoungsModulus = 1.0e6,
            PoissonRatio = 0.3,
            YieldStrength = 5.0e4,
            UltimateTensileStrength = 1.0e5,
            HardeningCoefficient = 0.0,
            Ductility = 0.8,
            ThermalConductivity = 0.03,
            SpecificHeat = 1500,
            MeltingPoint = 423
        };

        #endregion

        #region Game-Friendly Presets

        /// <summary>
        /// Creates a bouncy ball material for games.
        /// </summary>
        public static PhysicsMaterial BouncyBall() => new("BouncyBall")
        {
            Density = 1100,
            Restitution = 0.95,
            StaticFriction = 0.6,
            DynamicFriction = 0.5,
            YoungsModulus = 1.0e7,
            PoissonRatio = 0.49,
            YieldStrength = 1.0e8,
            UltimateTensileStrength = 2.0e8,
            Ductility = 0.95
        };

        /// <summary>
        /// Creates a frictionless material for ice-like surfaces in games.
        /// </summary>
        public static PhysicsMaterial Frictionless() => new("Frictionless")
        {
            Density = 1000,
            Restitution = 0.5,
            StaticFriction = 0.0,
            DynamicFriction = 0.0,
            YoungsModulus = 1.0e10,
            PoissonRatio = 0.3,
            YieldStrength = 1.0e10,
            UltimateTensileStrength = 1.0e10,
            Ductility = 0.0
        };

        /// <summary>
        /// Creates a super sticky material for games.
        /// </summary>
        public static PhysicsMaterial Sticky() => new("Sticky")
        {
            Density = 1200,
            Restitution = 0.0,
            StaticFriction = 2.0,
            DynamicFriction = 1.5,
            YoungsModulus = 1.0e6,
            PoissonRatio = 0.45,
            YieldStrength = 1.0e5,
            UltimateTensileStrength = 2.0e5,
            Ductility = 0.9
        };

        #endregion
    }
}
