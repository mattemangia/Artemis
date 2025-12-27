using System;
using System.Collections.Generic;
using Artemis.Core;
using Artemis.Bodies;
using Artemis.Materials;

namespace Artemis.Particles
{
    /// <summary>
    /// Manages combustion, fire spread, and fire-water interaction.
    /// </summary>
    public class CombustionSystem
    {
        #region Fields

        private readonly List<FireSimulation> _fires;
        private readonly List<FluidSimulation> _waterBodies;
        private readonly List<CombustibleObject> _combustibles;
        private readonly Random _random;

        #endregion

        #region Properties

        /// <summary>
        /// Gets or sets the fire spread rate (chance per second).
        /// </summary>
        public double FireSpreadRate { get; set; } = 0.5;

        /// <summary>
        /// Gets or sets the fire spread distance.
        /// </summary>
        public double FireSpreadDistance { get; set; } = 2.0;

        /// <summary>
        /// Gets or sets the water extinguish radius.
        /// </summary>
        public double WaterExtinguishRadius { get; set; } = 1.5;

        /// <summary>
        /// Gets or sets the minimum water density to extinguish fire.
        /// </summary>
        public double WaterExtinguishThreshold { get; set; } = 0.3;

        /// <summary>
        /// Gets all active fires.
        /// </summary>
        public IReadOnlyList<FireSimulation> Fires => _fires;

        /// <summary>
        /// Gets all water bodies.
        /// </summary>
        public IReadOnlyList<FluidSimulation> WaterBodies => _waterBodies;

        /// <summary>
        /// Gets all combustible objects.
        /// </summary>
        public IReadOnlyList<CombustibleObject> Combustibles => _combustibles;

        /// <summary>
        /// Event raised when an object catches fire.
        /// </summary>
        public event Action<CombustibleObject>? OnIgnition;

        /// <summary>
        /// Event raised when a fire is extinguished.
        /// </summary>
        public event Action<FireSimulation, ExtinguishReason>? OnExtinguish;

        /// <summary>
        /// Event raised when an object burns out completely.
        /// </summary>
        public event Action<CombustibleObject>? OnBurnedOut;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new combustion system.
        /// </summary>
        public CombustionSystem()
        {
            _fires = new List<FireSimulation>();
            _waterBodies = new List<FluidSimulation>();
            _combustibles = new List<CombustibleObject>();
            _random = new Random();
        }

        #endregion

        #region Registration

        /// <summary>
        /// Registers a fire with the combustion system.
        /// </summary>
        public void RegisterFire(FireSimulation fire)
        {
            if (!_fires.Contains(fire))
                _fires.Add(fire);
        }

        /// <summary>
        /// Unregisters a fire.
        /// </summary>
        public void UnregisterFire(FireSimulation fire)
        {
            _fires.Remove(fire);
        }

        /// <summary>
        /// Registers a water body (for extinguishing fires).
        /// </summary>
        public void RegisterWater(FluidSimulation water)
        {
            if (!_waterBodies.Contains(water))
                _waterBodies.Add(water);
        }

        /// <summary>
        /// Unregisters a water body.
        /// </summary>
        public void UnregisterWater(FluidSimulation water)
        {
            _waterBodies.Remove(water);
        }

        /// <summary>
        /// Registers a combustible object.
        /// </summary>
        public CombustibleObject RegisterCombustible(
            IPhysicsBody body,
            double flammability = 0.5,
            double burnTime = 10.0,
            double ignitionTemperature = 500.0)
        {
            var combustible = new CombustibleObject
            {
                Body = body,
                Flammability = flammability,
                MaxBurnTime = burnTime,
                RemainingBurnTime = burnTime,
                IgnitionTemperature = ignitionTemperature,
                IsOnFire = false
            };
            _combustibles.Add(combustible);
            return combustible;
        }

        /// <summary>
        /// Creates a combustible from a flammable material preset.
        /// </summary>
        public CombustibleObject RegisterCombustible(IPhysicsBody body, FlammableMaterialData material)
        {
            return RegisterCombustible(
                body,
                material.Flammability,
                material.BurnTime,
                material.IgnitionTemperature
            );
        }

        /// <summary>
        /// Unregisters a combustible object.
        /// </summary>
        public void UnregisterCombustible(CombustibleObject combustible)
        {
            _combustibles.Remove(combustible);
        }

        #endregion

        #region Update

        /// <summary>
        /// Updates the combustion system.
        /// </summary>
        /// <param name="deltaTime">Time step in seconds.</param>
        public void Update(double deltaTime)
        {
            // Update all fires
            foreach (var fire in _fires.ToArray()) // ToArray to allow removal during iteration
            {
                fire.Update(deltaTime);

                // Check if fire is extinguished naturally
                if (!fire.IsBurning && fire.Particles.AliveCount == 0)
                {
                    _fires.Remove(fire);
                    OnExtinguish?.Invoke(fire, ExtinguishReason.Natural);
                }
            }

            // Check water-fire interaction
            CheckWaterFireInteraction();

            // Check fire spread to combustibles
            CheckFireSpread(deltaTime);

            // Update burning objects
            UpdateBurningObjects(deltaTime);
        }

        private void CheckWaterFireInteraction()
        {
            foreach (var water in _waterBodies)
            {
                foreach (var fire in _fires.ToArray())
                {
                    if (!fire.IsBurning)
                        continue;

                    // Check if water particles are near fire
                    double waterNearFire = CountWaterNearPosition(
                        water,
                        fire.EmissionPosition,
                        WaterExtinguishRadius
                    );

                    if (waterNearFire > WaterExtinguishThreshold)
                    {
                        // Water extinguishes fire
                        fire.Intensity -= waterNearFire * 0.1;

                        if (fire.Intensity <= 0)
                        {
                            fire.Extinguish();
                            OnExtinguish?.Invoke(fire, ExtinguishReason.Water);

                            // Create steam effect
                            CreateSteamEffect(fire.EmissionPosition);
                        }
                    }

                    // Also check water hitting fire particles
                    ExtinguishFireParticlesWithWater(fire, water);
                }
            }
        }

        private double CountWaterNearPosition(FluidSimulation water, Vector3D position, double radius)
        {
            int count = 0;
            double radiusSq = radius * radius;

            foreach (var particle in water.Particles)
            {
                if (!particle.IsActive)
                    continue;

                if (Vector3D.DistanceSquared(particle.Position, position) < radiusSq)
                    count++;
            }

            return (double)count / 100.0; // Normalize
        }

        private void ExtinguishFireParticlesWithWater(FireSimulation fire, FluidSimulation water)
        {
            var fireParticles = fire.Particles.Particles;
            double checkRadius = 0.5;
            double checkRadiusSq = checkRadius * checkRadius;

            for (int i = 0; i < fire.Particles.MaxParticles; i++)
            {
                if (!fireParticles[i].IsAlive)
                    continue;

                // Check if any water particle is near this fire particle
                foreach (var waterParticle in water.Particles)
                {
                    if (!waterParticle.IsActive)
                        continue;

                    if (Vector3D.DistanceSquared(fireParticles[i].Position, waterParticle.Position) < checkRadiusSq)
                    {
                        // Water touches fire - kill the fire particle
                        fireParticles[i].Kill();

                        // Small steam puff
                        if (_random.NextDouble() < 0.2)
                        {
                            // Would spawn steam particle here
                        }
                        break;
                    }
                }
            }
        }

        private void CheckFireSpread(double deltaTime)
        {
            foreach (var fire in _fires)
            {
                if (!fire.IsBurning)
                    continue;

                foreach (var combustible in _combustibles)
                {
                    if (combustible.IsOnFire || combustible.IsBurnedOut)
                        continue;

                    // Check distance to fire
                    double distance = Vector3D.Distance(
                        fire.EmissionPosition,
                        combustible.Body.Position
                    );

                    if (distance < FireSpreadDistance)
                    {
                        // Chance to ignite based on distance and flammability
                        double igniteChance = FireSpreadRate *
                            combustible.Flammability *
                            (1.0 - distance / FireSpreadDistance) *
                            deltaTime;

                        if (_random.NextDouble() < igniteChance)
                        {
                            IgniteObject(combustible);
                        }
                    }
                }
            }
        }

        private void UpdateBurningObjects(double deltaTime)
        {
            foreach (var combustible in _combustibles.ToArray())
            {
                if (!combustible.IsOnFire)
                    continue;

                // Burn down
                combustible.RemainingBurnTime -= deltaTime;

                // Update fire intensity based on remaining burn time
                if (combustible.Fire != null)
                {
                    combustible.Fire.Intensity = combustible.RemainingBurnTime / combustible.MaxBurnTime;
                }

                // Check if burned out
                if (combustible.RemainingBurnTime <= 0)
                {
                    combustible.IsOnFire = false;
                    combustible.IsBurnedOut = true;

                    if (combustible.Fire != null)
                    {
                        combustible.Fire.Extinguish();
                        _fires.Remove(combustible.Fire);
                    }

                    OnBurnedOut?.Invoke(combustible);
                }
            }
        }

        #endregion

        #region Actions

        /// <summary>
        /// Ignites a combustible object.
        /// </summary>
        public FireSimulation? IgniteObject(CombustibleObject combustible)
        {
            if (combustible.IsOnFire || combustible.IsBurnedOut)
                return null;

            combustible.IsOnFire = true;

            // Create fire at object position
            var fire = new FireSimulation(combustible.Body.Position)
            {
                EmissionRadius = 0.3,
                Intensity = 1.0,
                FlameHeight = 1.0
            };

            combustible.Fire = fire;
            _fires.Add(fire);

            OnIgnition?.Invoke(combustible);
            return fire;
        }

        /// <summary>
        /// Ignites an object at a specific position.
        /// </summary>
        public FireSimulation StartFire(Vector3D position)
        {
            var fire = FireSimulation.Campfire(position);
            _fires.Add(fire);
            return fire;
        }

        /// <summary>
        /// Extinguishes all fires.
        /// </summary>
        public void ExtinguishAll()
        {
            foreach (var fire in _fires)
            {
                fire.Extinguish();
            }
            _fires.Clear();

            foreach (var combustible in _combustibles)
            {
                combustible.IsOnFire = false;
                combustible.Fire = null;
            }
        }

        /// <summary>
        /// Adds water at a position to extinguish fires.
        /// </summary>
        public void SplashWater(Vector3D position, double radius = 2.0)
        {
            foreach (var fire in _fires.ToArray())
            {
                double distance = Vector3D.Distance(fire.EmissionPosition, position);
                if (distance < radius)
                {
                    double effectiveness = 1.0 - distance / radius;
                    fire.Intensity -= effectiveness * 0.5;

                    if (fire.Intensity <= 0)
                    {
                        fire.Extinguish();
                        OnExtinguish?.Invoke(fire, ExtinguishReason.Water);
                        CreateSteamEffect(fire.EmissionPosition);
                    }
                }
            }
        }

        private void CreateSteamEffect(Vector3D position)
        {
            // Could create a steam simulation here
            // For now, this is a placeholder for the effect
        }

        #endregion

        #region Queries

        /// <summary>
        /// Gets the nearest fire to a position.
        /// </summary>
        public FireSimulation? GetNearestFire(Vector3D position)
        {
            FireSimulation? nearest = null;
            double nearestDist = double.MaxValue;

            foreach (var fire in _fires)
            {
                double dist = Vector3D.Distance(fire.EmissionPosition, position);
                if (dist < nearestDist)
                {
                    nearestDist = dist;
                    nearest = fire;
                }
            }

            return nearest;
        }

        /// <summary>
        /// Checks if a position is on fire.
        /// </summary>
        public bool IsPositionOnFire(Vector3D position, double radius = 1.0)
        {
            foreach (var fire in _fires)
            {
                if (Vector3D.Distance(fire.EmissionPosition, position) < radius)
                    return true;
            }
            return false;
        }

        /// <summary>
        /// Gets the temperature at a position.
        /// </summary>
        public double GetTemperatureAtPosition(Vector3D position)
        {
            double temp = 293; // Ambient ~20°C

            foreach (var fire in _fires)
            {
                if (!fire.IsBurning)
                    continue;

                double distance = Vector3D.Distance(fire.EmissionPosition, position);
                double fireRadius = fire.FlameHeight * 2;

                if (distance < fireRadius)
                {
                    // Temperature falls off with distance
                    double falloff = 1.0 - distance / fireRadius;
                    temp += fire.InitialTemperature * falloff * fire.Intensity;
                }
            }

            return temp;
        }

        #endregion
    }

    /// <summary>
    /// Represents a combustible object that can catch fire.
    /// </summary>
    public class CombustibleObject
    {
        /// <summary>
        /// The physics body.
        /// </summary>
        public IPhysicsBody Body { get; set; } = null!;

        /// <summary>
        /// Flammability (0-1). Higher = catches fire easier.
        /// </summary>
        public double Flammability { get; set; } = 0.5;

        /// <summary>
        /// Maximum burn time in seconds.
        /// </summary>
        public double MaxBurnTime { get; set; } = 10.0;

        /// <summary>
        /// Remaining burn time.
        /// </summary>
        public double RemainingBurnTime { get; set; }

        /// <summary>
        /// Temperature needed to ignite (Kelvin).
        /// </summary>
        public double IgnitionTemperature { get; set; } = 500.0;

        /// <summary>
        /// Whether currently on fire.
        /// </summary>
        public bool IsOnFire { get; set; }

        /// <summary>
        /// Whether completely burned out.
        /// </summary>
        public bool IsBurnedOut { get; set; }

        /// <summary>
        /// The fire attached to this object.
        /// </summary>
        public FireSimulation? Fire { get; set; }

        /// <summary>
        /// Gets the burn progress (0 = full, 1 = burned out).
        /// </summary>
        public double BurnProgress => 1.0 - RemainingBurnTime / MaxBurnTime;
    }

    /// <summary>
    /// Reason why a fire was extinguished.
    /// </summary>
    public enum ExtinguishReason
    {
        Natural,
        Water,
        Manual,
        BurnedOut
    }

    /// <summary>
    /// Preset flammable materials.
    /// </summary>
    public static class FlammableMaterial
    {
        public static readonly FlammableMaterialData Wood = new()
        {
            Name = "Wood",
            Flammability = 0.7,
            BurnTime = 30.0,
            IgnitionTemperature = 573 // ~300°C
        };

        public static readonly FlammableMaterialData Paper = new()
        {
            Name = "Paper",
            Flammability = 0.95,
            BurnTime = 5.0,
            IgnitionTemperature = 505 // ~232°C (Fahrenheit 451!)
        };

        public static readonly FlammableMaterialData Fabric = new()
        {
            Name = "Fabric",
            Flammability = 0.8,
            BurnTime = 10.0,
            IgnitionTemperature = 533
        };

        public static readonly FlammableMaterialData Oil = new()
        {
            Name = "Oil",
            Flammability = 0.9,
            BurnTime = 60.0,
            IgnitionTemperature = 473
        };

        public static readonly FlammableMaterialData Gasoline = new()
        {
            Name = "Gasoline",
            Flammability = 1.0,
            BurnTime = 20.0,
            IgnitionTemperature = 519 // ~246°C
        };

        public static readonly FlammableMaterialData Coal = new()
        {
            Name = "Coal",
            Flammability = 0.4,
            BurnTime = 120.0,
            IgnitionTemperature = 700
        };

        public static readonly FlammableMaterialData Leaves = new()
        {
            Name = "Leaves",
            Flammability = 0.85,
            BurnTime = 3.0,
            IgnitionTemperature = 450
        };

        public static readonly FlammableMaterialData Rubber = new()
        {
            Name = "Rubber",
            Flammability = 0.6,
            BurnTime = 45.0,
            IgnitionTemperature = 533
        };
    }

    /// <summary>
    /// Data for a flammable material.
    /// </summary>
    public struct FlammableMaterialData
    {
        public string Name;
        public double Flammability;
        public double BurnTime;
        public double IgnitionTemperature;
    }
}
