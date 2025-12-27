using System;
using System.Collections.Generic;
using System.Numerics;

namespace Artemis.Forces
{
    /// <summary>
    /// Types of propulsion systems
    /// </summary>
    public enum PropulsionType
    {
        Chemical,       // Traditional rocket engines (high thrust, low efficiency)
        Ion,            // Ion thrusters (low thrust, high efficiency)
        Nuclear,        // Nuclear thermal rockets
        Solid,          // Solid rocket boosters
        Monopropellant, // Single propellant thrusters (RCS)
        Bipropellant,   // Two-propellant engines
        Hybrid,         // Hybrid rocket engines
        Cold,           // Cold gas thrusters
        Hall,           // Hall-effect thrusters
        Photon          // Solar sails / photon drives
    }

    /// <summary>
    /// Represents fuel/propellant for a propulsion system
    /// </summary>
    public class Propellant
    {
        public string Name { get; set; } = "Generic";
        public float Mass { get; set; }           // kg
        public float MaxMass { get; set; }        // kg (tank capacity)
        public float Density { get; set; }        // kg/m³
        public float EnergyDensity { get; set; }  // J/kg
        public bool IsExhausted => Mass <= 0;
        public float FillRatio => MaxMass > 0 ? Mass / MaxMass : 0;

        public static Propellant LiquidOxygen(float mass) => new()
        {
            Name = "LOX",
            Mass = mass,
            MaxMass = mass,
            Density = 1141f,
            EnergyDensity = 13.4e6f
        };

        public static Propellant LiquidHydrogen(float mass) => new()
        {
            Name = "LH2",
            Mass = mass,
            MaxMass = mass,
            Density = 70.8f,
            EnergyDensity = 142e6f
        };

        public static Propellant RP1(float mass) => new()
        {
            Name = "RP-1",
            Mass = mass,
            MaxMass = mass,
            Density = 810f,
            EnergyDensity = 43e6f
        };

        public static Propellant Xenon(float mass) => new()
        {
            Name = "Xenon",
            Mass = mass,
            MaxMass = mass,
            Density = 2940f,
            EnergyDensity = 0  // Ion engines use electrical power
        };

        public static Propellant Hydrazine(float mass) => new()
        {
            Name = "N2H4",
            Mass = mass,
            MaxMass = mass,
            Density = 1021f,
            EnergyDensity = 1.6e6f
        };
    }

    /// <summary>
    /// Represents a rocket engine or thruster
    /// </summary>
    public class Engine
    {
        public string Name { get; set; } = "Engine";
        public PropulsionType Type { get; set; }
        public float ThrustVacuum { get; set; }       // N (Newtons)
        public float ThrustSeaLevel { get; set; }    // N (Newtons) - less due to atmosphere
        public float SpecificImpulse { get; set; }   // seconds (Isp)
        public float MassFlowRate { get; set; }      // kg/s
        public float DryMass { get; set; }           // kg
        public float GimbalRange { get; set; }       // degrees
        public float ThrottleMin { get; set; }       // 0-1 (minimum throttle, some engines can't throttle)
        public float ThrottleMax { get; set; } = 1f; // 0-1
        public float CurrentThrottle { get; set; }   // 0-1
        public bool IsActive { get; set; }
        public float RestartTime { get; set; }       // seconds to restart
        public int MaxRestarts { get; set; } = -1;   // -1 = unlimited
        public int RestartCount { get; set; }

        // Calculated from Isp: thrust = Isp * g0 * massFlowRate
        private const float G0 = 9.80665f;

        /// <summary>
        /// Get current thrust based on throttle and atmospheric pressure
        /// </summary>
        public float GetThrust(float atmosphericPressure = 0)
        {
            if (!IsActive || CurrentThrottle < ThrottleMin) return 0;

            // Interpolate between sea level and vacuum thrust based on pressure
            float thrustAtPressure = MathF.Max(0, ThrustVacuum - (ThrustVacuum - ThrustSeaLevel) * atmosphericPressure);
            return thrustAtPressure * Math.Clamp(CurrentThrottle, ThrottleMin, ThrottleMax);
        }

        /// <summary>
        /// Get current fuel consumption rate
        /// </summary>
        public float GetFuelConsumption()
        {
            if (!IsActive || CurrentThrottle < ThrottleMin) return 0;
            return MassFlowRate * Math.Clamp(CurrentThrottle, ThrottleMin, ThrottleMax);
        }

        /// <summary>
        /// Calculate exhaust velocity from Isp
        /// </summary>
        public float ExhaustVelocity => SpecificImpulse * G0;

        // Engine presets
        public static Engine Merlin1D() => new()
        {
            Name = "Merlin 1D",
            Type = PropulsionType.Bipropellant,
            ThrustVacuum = 981_000,
            ThrustSeaLevel = 845_000,
            SpecificImpulse = 311,
            MassFlowRate = 321,
            DryMass = 470,
            GimbalRange = 5,
            ThrottleMin = 0.39f
        };

        public static Engine Raptor() => new()
        {
            Name = "Raptor",
            Type = PropulsionType.Bipropellant,
            ThrustVacuum = 2_200_000,
            ThrustSeaLevel = 1_850_000,
            SpecificImpulse = 380,
            MassFlowRate = 590,
            DryMass = 1600,
            GimbalRange = 15,
            ThrottleMin = 0.4f
        };

        public static Engine RL10() => new()
        {
            Name = "RL-10",
            Type = PropulsionType.Bipropellant,
            ThrustVacuum = 110_000,
            ThrustSeaLevel = 0, // Vacuum only
            SpecificImpulse = 465,
            MassFlowRate = 24,
            DryMass = 167,
            GimbalRange = 4,
            ThrottleMin = 0.1f
        };

        public static Engine IonThruster() => new()
        {
            Name = "NEXT Ion",
            Type = PropulsionType.Ion,
            ThrustVacuum = 0.236f,
            ThrustSeaLevel = 0,
            SpecificImpulse = 4190,
            MassFlowRate = 5.76e-6f,
            DryMass = 12.7f,
            GimbalRange = 0,
            ThrottleMin = 0.1f
        };

        public static Engine HallThruster() => new()
        {
            Name = "SPT-140",
            Type = PropulsionType.Hall,
            ThrustVacuum = 0.29f,
            ThrustSeaLevel = 0,
            SpecificImpulse = 1770,
            MassFlowRate = 1.67e-5f,
            DryMass = 8.5f,
            GimbalRange = 0,
            ThrottleMin = 0.2f
        };

        public static Engine RCS() => new()
        {
            Name = "RCS Thruster",
            Type = PropulsionType.Monopropellant,
            ThrustVacuum = 440,
            ThrustSeaLevel = 400,
            SpecificImpulse = 220,
            MassFlowRate = 0.2f,
            DryMass = 2,
            GimbalRange = 0,
            ThrottleMin = 0
        };

        public static Engine SolidBooster() => new()
        {
            Name = "SRB",
            Type = PropulsionType.Solid,
            ThrustVacuum = 12_000_000,
            ThrustSeaLevel = 11_800_000,
            SpecificImpulse = 268,
            MassFlowRate = 4500,
            DryMass = 68000,
            GimbalRange = 0,
            ThrottleMin = 1f, // Cannot throttle
            MaxRestarts = 0   // Cannot restart
        };
    }

    /// <summary>
    /// Represents a complete spacecraft propulsion stage
    /// </summary>
    public class PropulsionStage
    {
        public string Name { get; set; } = "Stage";
        public List<Engine> Engines { get; } = new();
        public List<Propellant> Propellants { get; } = new();
        public float DryMass { get; set; }           // kg (structure, avionics, etc.)
        public float PayloadMass { get; set; }       // kg (cargo, next stage, etc.)
        public bool IsJettisoned { get; set; }
        public bool IsSeparated { get; set; }

        /// <summary>
        /// Total mass including propellants
        /// </summary>
        public float WetMass => DryMass + PayloadMass +
            Engines.Sum(e => e.DryMass) +
            Propellants.Sum(p => p.Mass);

        /// <summary>
        /// Mass without propellants
        /// </summary>
        public float EmptyMass => DryMass + PayloadMass + Engines.Sum(e => e.DryMass);

        /// <summary>
        /// Total propellant mass
        /// </summary>
        public float PropellantMass => Propellants.Sum(p => p.Mass);

        /// <summary>
        /// Mass ratio (wet/dry)
        /// </summary>
        public float MassRatio => EmptyMass > 0 ? WetMass / EmptyMass : 1;

        /// <summary>
        /// Total thrust from all active engines
        /// </summary>
        public float GetTotalThrust(float atmosphericPressure = 0)
        {
            return Engines.Where(e => e.IsActive).Sum(e => e.GetThrust(atmosphericPressure));
        }

        /// <summary>
        /// Calculate delta-v using Tsiolkovsky rocket equation
        /// </summary>
        public float CalculateDeltaV()
        {
            if (Engines.Count == 0) return 0;

            float avgIsp = Engines.Where(e => e.IsActive || Engines.All(x => !x.IsActive))
                .Average(e => e.SpecificImpulse);

            if (avgIsp <= 0) return 0;

            float exhaustVelocity = avgIsp * 9.80665f;
            return exhaustVelocity * MathF.Log(MassRatio);
        }

        /// <summary>
        /// Activate all engines
        /// </summary>
        public void ActivateEngines()
        {
            foreach (var engine in Engines)
            {
                if (engine.MaxRestarts < 0 || engine.RestartCount < engine.MaxRestarts)
                {
                    engine.IsActive = true;
                    engine.RestartCount++;
                }
            }
        }

        /// <summary>
        /// Shutdown all engines
        /// </summary>
        public void ShutdownEngines()
        {
            foreach (var engine in Engines)
            {
                engine.IsActive = false;
                engine.CurrentThrottle = 0;
            }
        }

        /// <summary>
        /// Set throttle for all engines
        /// </summary>
        public void SetThrottle(float throttle)
        {
            foreach (var engine in Engines)
            {
                engine.CurrentThrottle = throttle;
            }
        }
    }

    /// <summary>
    /// Complete spacecraft with multiple stages
    /// </summary>
    public class Spacecraft
    {
        public string Name { get; set; } = "Spacecraft";
        public List<PropulsionStage> Stages { get; } = new();
        public int CurrentStageIndex { get; set; }
        public Vector3 Position { get; set; }
        public Vector3 Velocity { get; set; }
        public Vector3 Orientation { get; set; }     // Euler angles
        public Vector3 AngularVelocity { get; set; }
        public float TotalDeltaVExpended { get; set; }

        public PropulsionStage? CurrentStage =>
            CurrentStageIndex >= 0 && CurrentStageIndex < Stages.Count
                ? Stages[CurrentStageIndex]
                : null;

        /// <summary>
        /// Total mass of spacecraft (all non-jettisoned stages)
        /// </summary>
        public float TotalMass => Stages.Where(s => !s.IsJettisoned).Sum(s => s.WetMass);

        /// <summary>
        /// Total delta-v remaining across all stages
        /// </summary>
        public float TotalDeltaV => Stages.Where(s => !s.IsJettisoned).Sum(s => s.CalculateDeltaV());

        /// <summary>
        /// Stage separation - jettison current stage and activate next
        /// </summary>
        public bool SeparateStage()
        {
            if (CurrentStage == null) return false;

            CurrentStage.IsSeparated = true;
            CurrentStage.IsJettisoned = true;
            CurrentStage.ShutdownEngines();

            CurrentStageIndex++;

            if (CurrentStage != null)
            {
                CurrentStage.ActivateEngines();
                return true;
            }
            return false;
        }

        /// <summary>
        /// Update spacecraft physics
        /// </summary>
        public void Update(float deltaTime, float atmosphericPressure = 0, Vector3? gravity = null)
        {
            if (CurrentStage == null) return;

            // Calculate thrust vector
            float thrust = CurrentStage.GetTotalThrust(atmosphericPressure);
            Vector3 thrustDirection = GetThrustDirection();
            Vector3 thrustForce = thrustDirection * thrust;

            // Consume propellant
            float fuelRate = CurrentStage.Engines.Where(e => e.IsActive).Sum(e => e.GetFuelConsumption());
            ConsumePropellant(fuelRate * deltaTime);

            // Calculate acceleration (F = ma)
            float mass = TotalMass;
            Vector3 acceleration = mass > 0 ? thrustForce / mass : Vector3.Zero;

            // Add gravity
            if (gravity.HasValue)
            {
                acceleration += gravity.Value;
            }

            // Integrate motion
            Velocity += acceleration * deltaTime;
            Position += Velocity * deltaTime;

            // Track delta-v expended
            float deltaV = (thrust / mass) * deltaTime;
            TotalDeltaVExpended += deltaV;

            // Auto-stage if out of fuel
            if (CurrentStage.PropellantMass <= 0 && CurrentStageIndex < Stages.Count - 1)
            {
                SeparateStage();
            }
        }

        private Vector3 GetThrustDirection()
        {
            // Convert Euler angles to direction vector
            float pitch = Orientation.X;
            float yaw = Orientation.Y;

            return new Vector3(
                MathF.Cos(pitch) * MathF.Sin(yaw),
                MathF.Sin(pitch),
                MathF.Cos(pitch) * MathF.Cos(yaw)
            );
        }

        private void ConsumePropellant(float amount)
        {
            if (CurrentStage == null) return;

            float remaining = amount;
            foreach (var prop in CurrentStage.Propellants)
            {
                if (prop.Mass > 0)
                {
                    float consumed = MathF.Min(prop.Mass, remaining);
                    prop.Mass -= consumed;
                    remaining -= consumed;
                    if (remaining <= 0) break;
                }
            }
        }

        // Preset spacecraft configurations
        public static Spacecraft Falcon9()
        {
            var sc = new Spacecraft { Name = "Falcon 9" };

            // First stage
            var stage1 = new PropulsionStage { Name = "Stage 1", DryMass = 22200 };
            for (int i = 0; i < 9; i++)
            {
                stage1.Engines.Add(Engine.Merlin1D());
            }
            stage1.Propellants.Add(Propellant.RP1(123500));
            stage1.Propellants.Add(Propellant.LiquidOxygen(287400));

            // Second stage
            var stage2 = new PropulsionStage { Name = "Stage 2", DryMass = 4000 };
            var mvac = Engine.Merlin1D();
            mvac.Name = "Merlin 1D Vacuum";
            mvac.SpecificImpulse = 348;
            mvac.ThrustVacuum = 934_000;
            mvac.ThrustSeaLevel = 0;
            stage2.Engines.Add(mvac);
            stage2.Propellants.Add(Propellant.RP1(32300));
            stage2.Propellants.Add(Propellant.LiquidOxygen(75200));

            sc.Stages.Add(stage1);
            sc.Stages.Add(stage2);
            sc.CurrentStageIndex = 0;
            stage1.ActivateEngines();

            return sc;
        }

        public static Spacecraft IonProbe()
        {
            var sc = new Spacecraft { Name = "Ion Probe" };

            var stage = new PropulsionStage { Name = "Primary", DryMass = 100 };
            stage.Engines.Add(Engine.IonThruster());
            stage.Propellants.Add(Propellant.Xenon(50));

            sc.Stages.Add(stage);
            sc.CurrentStageIndex = 0;
            stage.ActivateEngines();

            return sc;
        }

        public static Spacecraft SaturnV()
        {
            var sc = new Spacecraft { Name = "Saturn V" };

            // S-IC First Stage (5x F-1 engines)
            var stage1 = new PropulsionStage { Name = "S-IC", DryMass = 130000 };
            for (int i = 0; i < 5; i++)
            {
                var f1 = new Engine
                {
                    Name = "F-1",
                    Type = PropulsionType.Bipropellant,
                    ThrustVacuum = 7_740_000,
                    ThrustSeaLevel = 6_770_000,
                    SpecificImpulse = 304,
                    MassFlowRate = 2578,
                    DryMass = 8400,
                    GimbalRange = 6,
                    ThrottleMin = 0.65f
                };
                stage1.Engines.Add(f1);
            }
            stage1.Propellants.Add(Propellant.RP1(770000));
            stage1.Propellants.Add(Propellant.LiquidOxygen(1500000));

            // S-II Second Stage (5x J-2 engines)
            var stage2 = new PropulsionStage { Name = "S-II", DryMass = 36000 };
            for (int i = 0; i < 5; i++)
            {
                var j2 = new Engine
                {
                    Name = "J-2",
                    Type = PropulsionType.Bipropellant,
                    ThrustVacuum = 1_033_000,
                    ThrustSeaLevel = 0,
                    SpecificImpulse = 421,
                    MassFlowRate = 250,
                    DryMass = 1788,
                    GimbalRange = 7.5f,
                    ThrottleMin = 0.77f
                };
                stage2.Engines.Add(j2);
            }
            stage2.Propellants.Add(Propellant.LiquidHydrogen(72000));
            stage2.Propellants.Add(Propellant.LiquidOxygen(360000));

            // S-IVB Third Stage (1x J-2 engine)
            var stage3 = new PropulsionStage { Name = "S-IVB", DryMass = 11500 };
            stage3.Engines.Add(new Engine
            {
                Name = "J-2",
                Type = PropulsionType.Bipropellant,
                ThrustVacuum = 1_033_000,
                ThrustSeaLevel = 0,
                SpecificImpulse = 421,
                MassFlowRate = 250,
                DryMass = 1788,
                GimbalRange = 7.5f,
                ThrottleMin = 0.77f,
                MaxRestarts = 2
            });
            stage3.Propellants.Add(Propellant.LiquidHydrogen(19000));
            stage3.Propellants.Add(Propellant.LiquidOxygen(87000));

            sc.Stages.Add(stage1);
            sc.Stages.Add(stage2);
            sc.Stages.Add(stage3);
            sc.CurrentStageIndex = 0;
            stage1.ActivateEngines();

            return sc;
        }
    }

    /// <summary>
    /// Manages spacecraft propulsion physics integration
    /// </summary>
    public class PropulsionSystem
    {
        private readonly List<Spacecraft> _spacecraft = new();
        private readonly OrbitalMechanics _orbitalMechanics = new();

        public IReadOnlyList<Spacecraft> Spacecraft => _spacecraft;

        /// <summary>
        /// Add a spacecraft to the simulation
        /// </summary>
        public void AddSpacecraft(Spacecraft spacecraft)
        {
            _spacecraft.Add(spacecraft);
        }

        /// <summary>
        /// Remove a spacecraft from the simulation
        /// </summary>
        public void RemoveSpacecraft(Spacecraft spacecraft)
        {
            _spacecraft.Remove(spacecraft);
        }

        /// <summary>
        /// Update all spacecraft
        /// </summary>
        public void Update(float deltaTime, List<CelestialBody>? celestialBodies = null)
        {
            foreach (var sc in _spacecraft)
            {
                // Calculate gravity from celestial bodies
                Vector3 gravity = Vector3.Zero;
                if (celestialBodies != null)
                {
                    foreach (var body in celestialBodies)
                    {
                        Vector3 toBody = (Vector3)body.Position - sc.Position;
                        float distance = toBody.Length();
                        if (distance > 0)
                        {
                        float gravMag = (float)(body.Mu / (distance * distance));
                            gravity += Vector3.Normalize(toBody) * gravMag;
                        }
                    }
                }

                // Calculate atmospheric pressure (simple exponential model)
                float atmosphericPressure = 0;
                if (celestialBodies != null)
                {
                    foreach (var body in celestialBodies)
                    {
                        float altitude = (sc.Position - (Vector3)body.Position).Length() - (float)body.Radius;
                        if (altitude < (float)body.AtmosphereHeight && altitude >= 0)
                        {
                            // Exponential atmosphere model
                            float scaleHeight = (float)body.AtmosphereHeight / 7f; // Approximate scale height
                            atmosphericPressure = MathF.Max(atmosphericPressure,
                                MathF.Exp(-altitude / scaleHeight));
                        }
                    }
                }

                sc.Update(deltaTime, atmosphericPressure, gravity);
            }
        }

        /// <summary>
        /// Calculate burn time for a maneuver
        /// </summary>
        public static float CalculateBurnTime(Spacecraft spacecraft, float deltaV)
        {
            var stage = spacecraft.CurrentStage;
            if (stage == null) return float.PositiveInfinity;

            float thrust = stage.GetTotalThrust(0);
            float mass = spacecraft.TotalMass;
            float massFlowRate = stage.Engines.Where(e => e.IsActive).Sum(e => e.MassFlowRate);

            if (thrust <= 0 || massFlowRate <= 0) return float.PositiveInfinity;

            // Use Tsiolkovsky equation to find final mass, then calculate time
            float avgIsp = stage.Engines.Where(e => e.IsActive).Average(e => e.SpecificImpulse);
            float exhaustVelocity = avgIsp * 9.80665f;
            float finalMass = mass / MathF.Exp(deltaV / exhaustVelocity);
            float propellantNeeded = mass - finalMass;

            return propellantNeeded / massFlowRate;
        }

        /// <summary>
        /// Execute a maneuver with automatic staging
        /// </summary>
        public void ExecuteManeuver(Spacecraft spacecraft, Vector3 deltaVVector, float deltaTime)
        {
            float targetDeltaV = deltaVVector.Length();
            float deltaVRemaining = targetDeltaV;
            Vector3 direction = Vector3.Normalize(deltaVVector);

            // Set orientation to burn direction
            spacecraft.Orientation = new Vector3(
                MathF.Asin(direction.Y),
                MathF.Atan2(direction.X, direction.Z),
                0
            );

            while (deltaVRemaining > 0 && spacecraft.CurrentStage != null)
            {
                var stage = spacecraft.CurrentStage;
                stage.SetThrottle(1f);

                float thrust = stage.GetTotalThrust(0);
                float mass = spacecraft.TotalMass;

                if (thrust <= 0)
                {
                    spacecraft.SeparateStage();
                    continue;
                }

                float acceleration = thrust / mass;
                float deltaVThisStep = acceleration * deltaTime;

                if (deltaVThisStep >= deltaVRemaining)
                {
                    // Final burn - throttle down
                    float burnFraction = deltaVRemaining / deltaVThisStep;
                    stage.SetThrottle(burnFraction);
                    deltaVRemaining = 0;
                }
                else
                {
                    deltaVRemaining -= deltaVThisStep;
                }

                spacecraft.Update(deltaTime, 0, null);

                // Check if stage is exhausted
                if (stage.PropellantMass <= 0)
                {
                    spacecraft.SeparateStage();
                }
            }

            // Shutdown engines after maneuver
            spacecraft.CurrentStage?.ShutdownEngines();
        }
    }
}
