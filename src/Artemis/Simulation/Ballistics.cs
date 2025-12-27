using System;
using System.Collections.Generic;
using System.Numerics;
using Artemis.Core;

namespace Artemis.Simulation
{
    /// <summary>
    /// Bullet/projectile configuration.
    /// </summary>
    public class Projectile
    {
        /// <summary>Unique identifier.</summary>
        public string Id { get; set; } = Guid.NewGuid().ToString();

        /// <summary>Projectile name.</summary>
        public string Name { get; set; } = "Projectile";

        /// <summary>Mass in kilograms.</summary>
        public float Mass { get; set; } = 0.01f;

        /// <summary>Diameter in meters.</summary>
        public float Diameter { get; set; } = 0.009f;

        /// <summary>Drag coefficient (typical bullet: 0.3-0.5).</summary>
        public float DragCoefficient { get; set; } = 0.4f;

        /// <summary>Ballistic coefficient (higher = less drag, more stable).</summary>
        public float BallisticCoefficient => Mass / (DragCoefficient * CrossSectionalArea);

        /// <summary>Cross-sectional area in m².</summary>
        public float CrossSectionalArea => MathF.PI * (Diameter / 2) * (Diameter / 2);

        /// <summary>Muzzle velocity in m/s.</summary>
        public float MuzzleVelocity { get; set; } = 800f;

        /// <summary>Muzzle energy in Joules.</summary>
        public float MuzzleEnergy => 0.5f * Mass * MuzzleVelocity * MuzzleVelocity;

        /// <summary>Spin rate in rad/s (for gyroscopic stability).</summary>
        public float SpinRate { get; set; } = 0f;

        /// <summary>Rifling twist rate (inches per turn, 0 = smoothbore).</summary>
        public float TwistRate { get; set; } = 0f;

        #region Presets

        /// <summary>9mm Parabellum (pistol).</summary>
        public static Projectile Bullet9mm() => new()
        {
            Name = "9mm Parabellum",
            Mass = 0.008f,          // 8g (124 grain)
            Diameter = 0.009f,      // 9mm
            DragCoefficient = 0.365f,
            MuzzleVelocity = 360f,  // ~1180 fps
            TwistRate = 10f
        };

        /// <summary>5.56x45mm NATO (rifle).</summary>
        public static Projectile Bullet556() => new()
        {
            Name = "5.56x45mm NATO",
            Mass = 0.004f,          // 4g (62 grain)
            Diameter = 0.00566f,    // 5.66mm
            DragCoefficient = 0.295f,
            MuzzleVelocity = 940f,  // ~3080 fps
            TwistRate = 7f
        };

        /// <summary>7.62x51mm NATO (rifle).</summary>
        public static Projectile Bullet762() => new()
        {
            Name = "7.62x51mm NATO",
            Mass = 0.0095f,         // 9.5g (147 grain)
            Diameter = 0.00762f,    // 7.62mm
            DragCoefficient = 0.393f,
            MuzzleVelocity = 850f,  // ~2800 fps
            TwistRate = 12f
        };

        /// <summary>.50 BMG (heavy machine gun).</summary>
        public static Projectile Bullet50BMG() => new()
        {
            Name = ".50 BMG",
            Mass = 0.0427f,         // 42.7g (660 grain)
            Diameter = 0.0127f,     // 12.7mm
            DragCoefficient = 0.62f,
            MuzzleVelocity = 928f,  // ~3040 fps
            TwistRate = 15f
        };

        /// <summary>12 gauge shotgun slug.</summary>
        public static Projectile ShotgunSlug() => new()
        {
            Name = "12ga Slug",
            Mass = 0.028f,          // 28g (1 oz)
            Diameter = 0.0185f,     // 18.5mm
            DragCoefficient = 0.8f,
            MuzzleVelocity = 450f,  // ~1500 fps
            TwistRate = 0f          // Smoothbore
        };

        /// <summary>Arrow (compound bow).</summary>
        public static Projectile Arrow() => new()
        {
            Name = "Arrow",
            Mass = 0.025f,          // 25g
            Diameter = 0.008f,      // 8mm shaft
            DragCoefficient = 0.5f,
            MuzzleVelocity = 100f,  // ~330 fps
            SpinRate = 0f
        };

        /// <summary>Crossbow bolt.</summary>
        public static Projectile CrossbowBolt() => new()
        {
            Name = "Crossbow Bolt",
            Mass = 0.022f,          // 22g
            Diameter = 0.009f,
            DragCoefficient = 0.45f,
            MuzzleVelocity = 130f,  // ~425 fps
            SpinRate = 0f
        };

        /// <summary>Tank shell (APFSDS).</summary>
        public static Projectile TankShell() => new()
        {
            Name = "120mm APFSDS",
            Mass = 4.6f,            // 4.6kg penetrator
            Diameter = 0.027f,      // 27mm dart
            DragCoefficient = 0.15f,
            MuzzleVelocity = 1750f, // ~5740 fps
            TwistRate = 0f          // Fin-stabilized
        };

        /// <summary>Artillery shell (155mm).</summary>
        public static Projectile ArtilleryShell() => new()
        {
            Name = "155mm HE",
            Mass = 43.5f,           // 43.5kg
            Diameter = 0.155f,      // 155mm
            DragCoefficient = 0.3f,
            MuzzleVelocity = 827f,  // ~2710 fps
            SpinRate = 15000f       // RPM from rifling
        };

        /// <summary>Baseball (for comparison).</summary>
        public static Projectile Baseball() => new()
        {
            Name = "Baseball",
            Mass = 0.145f,          // 145g
            Diameter = 0.074f,      // 74mm
            DragCoefficient = 0.35f,
            MuzzleVelocity = 45f,   // ~100 mph pitch
            SpinRate = 200f         // Typical spin
        };

        /// <summary>Golf ball.</summary>
        public static Projectile GolfBall() => new()
        {
            Name = "Golf Ball",
            Mass = 0.0459f,         // 45.9g
            Diameter = 0.0427f,     // 42.7mm
            DragCoefficient = 0.25f, // Dimples reduce drag
            MuzzleVelocity = 70f,   // ~156 mph driver
            SpinRate = 300f         // Backspin
        };

        #endregion
    }

    /// <summary>
    /// Environment conditions affecting ballistics.
    /// </summary>
    public class BallisticEnvironment
    {
        /// <summary>Air density in kg/m³ (sea level: 1.225).</summary>
        public float AirDensity { get; set; } = 1.225f;

        /// <summary>Temperature in Celsius.</summary>
        public float Temperature { get; set; } = 15f;

        /// <summary>Pressure in hPa.</summary>
        public float Pressure { get; set; } = 1013.25f;

        /// <summary>Humidity 0-1.</summary>
        public float Humidity { get; set; } = 0.5f;

        /// <summary>Altitude in meters.</summary>
        public float Altitude { get; set; } = 0f;

        /// <summary>Wind velocity (x=right, y=up, z=forward relative to shooter).</summary>
        public Vector3 Wind { get; set; } = Vector3.Zero;

        /// <summary>Gravity vector.</summary>
        public Vector3 Gravity { get; set; } = new Vector3(0, -9.81f, 0);

        /// <summary>Speed of sound in m/s.</summary>
        public float SpeedOfSound => 331.3f * MathF.Sqrt(1 + Temperature / 273.15f);

        /// <summary>
        /// Calculates air density from temperature, pressure, and humidity.
        /// </summary>
        public void RecalculateAirDensity()
        {
            // Ideal gas law with humidity correction
            float tempK = Temperature + 273.15f;
            float satVaporPressure = 6.1078f * MathF.Pow(10, 7.5f * Temperature / (Temperature + 237.3f));
            float vaporPressure = Humidity * satVaporPressure;
            float dryPressure = Pressure - vaporPressure;

            // Rd = 287.05, Rv = 461.495
            AirDensity = (dryPressure * 100 / (287.05f * tempK)) +
                        (vaporPressure * 100 / (461.495f * tempK));
        }

        /// <summary>Standard sea level conditions.</summary>
        public static BallisticEnvironment Standard() => new();

        /// <summary>High altitude conditions (3000m).</summary>
        public static BallisticEnvironment HighAltitude() => new()
        {
            Altitude = 3000f,
            AirDensity = 0.9f,
            Temperature = 0f,
            Pressure = 700f
        };

        /// <summary>Hot desert conditions.</summary>
        public static BallisticEnvironment Desert() => new()
        {
            Temperature = 45f,
            Humidity = 0.1f,
            AirDensity = 1.1f
        };

        /// <summary>Arctic conditions.</summary>
        public static BallisticEnvironment Arctic() => new()
        {
            Temperature = -30f,
            Humidity = 0.3f,
            AirDensity = 1.45f
        };
    }

    /// <summary>
    /// Single point in a trajectory.
    /// </summary>
    public struct TrajectoryPoint
    {
        /// <summary>Time since launch in seconds.</summary>
        public float Time;

        /// <summary>Position in world space.</summary>
        public Vector3 Position;

        /// <summary>Velocity in world space.</summary>
        public Vector3 Velocity;

        /// <summary>Speed (magnitude of velocity).</summary>
        public float Speed => Velocity.Length();

        /// <summary>Kinetic energy in Joules.</summary>
        public float KineticEnergy;

        /// <summary>Distance traveled from origin.</summary>
        public float Distance;

        /// <summary>Mach number (velocity / speed of sound).</summary>
        public float MachNumber;

        /// <summary>Spin rate (decays over time).</summary>
        public float SpinRate;
    }

    /// <summary>
    /// Complete trajectory result.
    /// </summary>
    public class TrajectoryResult
    {
        /// <summary>All trajectory points.</summary>
        public List<TrajectoryPoint> Points { get; } = new();

        /// <summary>Total flight time.</summary>
        public float TotalTime { get; set; }

        /// <summary>Total distance traveled.</summary>
        public float TotalDistance { get; set; }

        /// <summary>Maximum height reached.</summary>
        public float MaxHeight { get; set; }

        /// <summary>Final position.</summary>
        public Vector3 FinalPosition { get; set; }

        /// <summary>Final velocity.</summary>
        public Vector3 FinalVelocity { get; set; }

        /// <summary>Final kinetic energy.</summary>
        public float FinalEnergy { get; set; }

        /// <summary>Impact angle in degrees.</summary>
        public float ImpactAngle { get; set; }

        /// <summary>Whether projectile went supersonic to subsonic.</summary>
        public bool WentSubsonic { get; set; }

        /// <summary>Time when projectile went subsonic.</summary>
        public float SubsonicTime { get; set; }

        /// <summary>Wind drift (lateral displacement).</summary>
        public float WindDrift { get; set; }

        /// <summary>Bullet drop from line of sight.</summary>
        public float Drop { get; set; }

        /// <summary>
        /// Gets position at specific time via interpolation.
        /// </summary>
        public Vector3 GetPositionAtTime(float time)
        {
            if (Points.Count == 0) return Vector3.Zero;
            if (time <= 0) return Points[0].Position;
            if (time >= TotalTime) return FinalPosition;

            for (int i = 1; i < Points.Count; i++)
            {
                if (Points[i].Time >= time)
                {
                    float t = (time - Points[i - 1].Time) / (Points[i].Time - Points[i - 1].Time);
                    return Vector3.Lerp(Points[i - 1].Position, Points[i].Position, t);
                }
            }

            return FinalPosition;
        }

        /// <summary>
        /// Gets time of flight to reach a distance.
        /// </summary>
        public float GetTimeToDistance(float distance)
        {
            for (int i = 1; i < Points.Count; i++)
            {
                if (Points[i].Distance >= distance)
                {
                    float t = (distance - Points[i - 1].Distance) /
                             (Points[i].Distance - Points[i - 1].Distance);
                    return Points[i - 1].Time + t * (Points[i].Time - Points[i - 1].Time);
                }
            }
            return TotalTime;
        }
    }

    /// <summary>
    /// Ballistics calculation system for projectile trajectory simulation.
    /// </summary>
    public class BallisticsSystem
    {
        /// <summary>Default time step for trajectory calculation.</summary>
        public float TimeStep { get; set; } = 0.001f;

        /// <summary>Maximum simulation time in seconds.</summary>
        public float MaxTime { get; set; } = 10f;

        /// <summary>Minimum velocity to continue simulation.</summary>
        public float MinVelocity { get; set; } = 10f;

        /// <summary>Ground plane Y level (stops simulation on impact).</summary>
        public float GroundLevel { get; set; } = 0f;

        /// <summary>Whether to simulate spin decay.</summary>
        public bool SimulateSpinDecay { get; set; } = true;

        /// <summary>Whether to use advanced drag model (Mach-dependent).</summary>
        public bool UseAdvancedDragModel { get; set; } = true;

        /// <summary>
        /// Calculates the full trajectory of a projectile.
        /// </summary>
        /// <param name="projectile">Projectile configuration.</param>
        /// <param name="origin">Starting position.</param>
        /// <param name="direction">Launch direction (will be normalized).</param>
        /// <param name="environment">Environment conditions.</param>
        /// <returns>Complete trajectory result.</returns>
        public TrajectoryResult CalculateTrajectory(
            Projectile projectile,
            Vector3 origin,
            Vector3 direction,
            BallisticEnvironment? environment = null)
        {
            environment ??= BallisticEnvironment.Standard();
            direction = Vector3.Normalize(direction);

            var result = new TrajectoryResult();
            var position = origin;
            var velocity = direction * projectile.MuzzleVelocity;
            var spinRate = projectile.SpinRate > 0 ? projectile.SpinRate :
                (projectile.TwistRate > 0 ? CalculateSpinRate(projectile) : 0);

            float time = 0;
            float distance = 0;
            float maxHeight = origin.Y;
            bool wasSupersonic = velocity.Length() > environment.SpeedOfSound;
            float speedOfSound = environment.SpeedOfSound;

            // Record initial point
            result.Points.Add(new TrajectoryPoint
            {
                Time = 0,
                Position = position,
                Velocity = velocity,
                KineticEnergy = 0.5f * projectile.Mass * velocity.LengthSquared(),
                Distance = 0,
                MachNumber = velocity.Length() / speedOfSound,
                SpinRate = spinRate
            });

            while (time < MaxTime)
            {
                float speed = velocity.Length();
                if (speed < MinVelocity) break;
                if (position.Y < GroundLevel) break;

                // Calculate drag force
                float dragCoeff = UseAdvancedDragModel
                    ? GetMachAdjustedDragCoefficient(projectile.DragCoefficient, speed / speedOfSound)
                    : projectile.DragCoefficient;

                // Drag force: F = 0.5 * rho * v² * Cd * A
                float dragMagnitude = 0.5f * environment.AirDensity * speed * speed *
                                     dragCoeff * projectile.CrossSectionalArea;
                Vector3 dragForce = -Vector3.Normalize(velocity) * dragMagnitude;

                // Wind effect (relative velocity)
                Vector3 relativeWind = velocity - environment.Wind;
                float windDragMagnitude = 0.5f * environment.AirDensity *
                                         relativeWind.LengthSquared() * dragCoeff *
                                         projectile.CrossSectionalArea;
                Vector3 windEffect = -Vector3.Normalize(relativeWind) * windDragMagnitude - dragForce;

                // Total acceleration
                Vector3 acceleration = environment.Gravity +
                                       (dragForce + windEffect) / projectile.Mass;

                // Spin decay
                if (SimulateSpinDecay && spinRate > 0)
                {
                    spinRate *= 1 - (0.01f * TimeStep); // Gradual decay
                }

                // Integrate (Euler method - could use RK4 for more accuracy)
                velocity += acceleration * TimeStep;
                position += velocity * TimeStep;
                time += TimeStep;

                float stepDistance = velocity.Length() * TimeStep;
                distance += stepDistance;

                maxHeight = MathF.Max(maxHeight, position.Y);

                // Check supersonic transition
                if (wasSupersonic && speed < speedOfSound)
                {
                    result.WentSubsonic = true;
                    result.SubsonicTime = time;
                    wasSupersonic = false;
                }

                // Record point (at reduced frequency for memory efficiency)
                if (result.Points.Count < 10000 && (result.Points.Count == 0 ||
                    time - result.Points[^1].Time >= TimeStep * 10))
                {
                    result.Points.Add(new TrajectoryPoint
                    {
                        Time = time,
                        Position = position,
                        Velocity = velocity,
                        KineticEnergy = 0.5f * projectile.Mass * speed * speed,
                        Distance = distance,
                        MachNumber = speed / speedOfSound,
                        SpinRate = spinRate
                    });
                }
            }

            // Finalize result
            result.TotalTime = time;
            result.TotalDistance = distance;
            result.MaxHeight = maxHeight;
            result.FinalPosition = position;
            result.FinalVelocity = velocity;
            result.FinalEnergy = 0.5f * projectile.Mass * velocity.LengthSquared();
            result.ImpactAngle = MathF.Atan2(-velocity.Y, MathF.Sqrt(velocity.X * velocity.X + velocity.Z * velocity.Z)) * 180 / MathF.PI;
            result.WindDrift = position.X - origin.X;
            result.Drop = origin.Y - position.Y;

            return result;
        }

        /// <summary>
        /// Calculates trajectory for hitting a specific target.
        /// </summary>
        /// <param name="projectile">Projectile configuration.</param>
        /// <param name="origin">Starting position.</param>
        /// <param name="target">Target position.</param>
        /// <param name="environment">Environment conditions.</param>
        /// <returns>Required launch angle in degrees, or null if impossible.</returns>
        public (float elevation, float windage)? CalculateFiringSolution(
            Projectile projectile,
            Vector3 origin,
            Vector3 target,
            BallisticEnvironment? environment = null)
        {
            environment ??= BallisticEnvironment.Standard();

            // Calculate flat distance and height difference
            Vector3 toTarget = target - origin;
            float horizontalDist = MathF.Sqrt(toTarget.X * toTarget.X + toTarget.Z * toTarget.Z);
            float heightDiff = toTarget.Y;

            // Initial guess using simple ballistic formula
            float v = projectile.MuzzleVelocity;
            float g = -environment.Gravity.Y;

            // Try to solve analytically for initial guess
            float angle1, angle2;
            if (!SolveBallisticAngle(v, g, horizontalDist, heightDiff, out angle1, out angle2))
            {
                return null; // Target unreachable
            }

            // Use lower angle (direct fire) as starting point
            float elevation = angle1;

            // Iterative refinement for wind and drag
            for (int i = 0; i < 20; i++)
            {
                // Calculate horizontal direction
                float azimuth = MathF.Atan2(toTarget.X, toTarget.Z);

                // Create direction vector
                Vector3 direction = new Vector3(
                    MathF.Sin(azimuth) * MathF.Cos(elevation),
                    MathF.Sin(elevation),
                    MathF.Cos(azimuth) * MathF.Cos(elevation)
                );

                // Simulate trajectory
                var result = CalculateTrajectory(projectile, origin, direction, environment);

                // Check impact point
                float impactDist = (result.FinalPosition - target).Length();
                if (impactDist < 0.1f) // Within 10cm
                {
                    float windage = MathF.Atan2(result.WindDrift, horizontalDist) * 180 / MathF.PI;
                    return (elevation * 180 / MathF.PI, windage);
                }

                // Adjust angle based on error
                float vertError = result.FinalPosition.Y - target.Y;
                float horzError = Vector3.Distance(
                    new Vector3(result.FinalPosition.X, 0, result.FinalPosition.Z),
                    new Vector3(target.X, 0, target.Z)) - horizontalDist;

                elevation += vertError * 0.001f;
            }

            return (elevation * 180 / MathF.PI, 0);
        }

        /// <summary>
        /// Predicts where a projectile will land.
        /// </summary>
        public Vector3 PredictImpact(
            Projectile projectile,
            Vector3 origin,
            Vector3 direction,
            BallisticEnvironment? environment = null)
        {
            var result = CalculateTrajectory(projectile, origin, direction, environment);
            return result.FinalPosition;
        }

        /// <summary>
        /// Calculates time of flight to a target distance.
        /// </summary>
        public float CalculateTimeOfFlight(
            Projectile projectile,
            float distance,
            float launchAngle,
            BallisticEnvironment? environment = null)
        {
            var direction = new Vector3(0, MathF.Sin(launchAngle * MathF.PI / 180),
                                        MathF.Cos(launchAngle * MathF.PI / 180));
            var result = CalculateTrajectory(projectile, Vector3.Zero, direction, environment);
            return result.GetTimeToDistance(distance);
        }

        /// <summary>
        /// Calculates bullet drop at a given distance.
        /// </summary>
        public float CalculateDrop(
            Projectile projectile,
            float distance,
            float launchAngle = 0,
            BallisticEnvironment? environment = null)
        {
            var direction = new Vector3(0, MathF.Sin(launchAngle * MathF.PI / 180),
                                        MathF.Cos(launchAngle * MathF.PI / 180));
            var result = CalculateTrajectory(projectile, Vector3.Zero, direction, environment);

            // Find position at target distance
            var pos = result.GetPositionAtTime(result.GetTimeToDistance(distance));
            return -pos.Y; // Positive drop = below line of sight
        }

        /// <summary>
        /// Calculates kinetic energy at a given distance.
        /// </summary>
        public float CalculateEnergyAtDistance(
            Projectile projectile,
            float distance,
            BallisticEnvironment? environment = null)
        {
            var result = CalculateTrajectory(projectile, Vector3.Zero, Vector3.UnitZ, environment);
            float time = result.GetTimeToDistance(distance);

            for (int i = 1; i < result.Points.Count; i++)
            {
                if (result.Points[i].Time >= time)
                {
                    float t = (time - result.Points[i - 1].Time) /
                             (result.Points[i].Time - result.Points[i - 1].Time);
                    return result.Points[i - 1].KineticEnergy +
                           t * (result.Points[i].KineticEnergy - result.Points[i - 1].KineticEnergy);
                }
            }

            return result.FinalEnergy;
        }

        /// <summary>
        /// Calculates maximum range for a projectile.
        /// </summary>
        public float CalculateMaxRange(
            Projectile projectile,
            BallisticEnvironment? environment = null)
        {
            environment ??= BallisticEnvironment.Standard();

            float maxRange = 0;
            float optimalAngle = 45; // Start with vacuum optimal

            // Search for optimal angle
            for (float angle = 30; angle <= 60; angle += 1)
            {
                var direction = new Vector3(0, MathF.Sin(angle * MathF.PI / 180),
                                            MathF.Cos(angle * MathF.PI / 180));
                var result = CalculateTrajectory(projectile, Vector3.Zero, direction, environment);

                float range = MathF.Sqrt(result.FinalPosition.X * result.FinalPosition.X +
                                        result.FinalPosition.Z * result.FinalPosition.Z);

                if (range > maxRange)
                {
                    maxRange = range;
                    optimalAngle = angle;
                }
            }

            return maxRange;
        }

        /// <summary>
        /// Creates a range table for a projectile.
        /// </summary>
        public Dictionary<float, (float drop, float time, float energy, float velocity)> CreateRangeTable(
            Projectile projectile,
            float maxDistance,
            float interval = 100f,
            BallisticEnvironment? environment = null)
        {
            var table = new Dictionary<float, (float, float, float, float)>();
            var result = CalculateTrajectory(projectile, Vector3.Zero, Vector3.UnitZ, environment);

            for (float dist = interval; dist <= maxDistance; dist += interval)
            {
                float time = result.GetTimeToDistance(dist);
                if (time >= result.TotalTime) break;

                var pos = result.GetPositionAtTime(time);
                float drop = -pos.Y;

                // Find energy and velocity at this time
                float energy = 0, velocity = 0;
                for (int i = 1; i < result.Points.Count; i++)
                {
                    if (result.Points[i].Time >= time)
                    {
                        float t = (time - result.Points[i - 1].Time) /
                                 (result.Points[i].Time - result.Points[i - 1].Time);
                        energy = result.Points[i - 1].KineticEnergy +
                                t * (result.Points[i].KineticEnergy - result.Points[i - 1].KineticEnergy);
                        velocity = result.Points[i - 1].Speed +
                                  t * (result.Points[i].Speed - result.Points[i - 1].Speed);
                        break;
                    }
                }

                table[dist] = (drop, time, energy, velocity);
            }

            return table;
        }

        #region Private Methods

        private float CalculateSpinRate(Projectile projectile)
        {
            if (projectile.TwistRate <= 0) return 0;

            // Convert twist rate (inches per turn) to rad/s
            // Spin = velocity / (twist_rate_inches * 0.0254) * 2π
            float twistMeters = projectile.TwistRate * 0.0254f;
            float turnsPerSecond = projectile.MuzzleVelocity / twistMeters;
            return turnsPerSecond * 2 * MathF.PI;
        }

        private float GetMachAdjustedDragCoefficient(float baseCd, float mach)
        {
            // Simplified drag coefficient model based on Mach number
            // Real models use G1, G7 drag functions with lookup tables

            if (mach < 0.8f)
            {
                // Subsonic - relatively constant
                return baseCd;
            }
            else if (mach < 1.2f)
            {
                // Transonic - drag increases dramatically
                float t = (mach - 0.8f) / 0.4f;
                return baseCd * (1 + 0.5f * t);
            }
            else if (mach < 2.0f)
            {
                // Low supersonic - drag decreases
                float t = (mach - 1.2f) / 0.8f;
                return baseCd * (1.5f - 0.3f * t);
            }
            else
            {
                // High supersonic
                return baseCd * 1.2f / MathF.Sqrt(mach);
            }
        }

        private bool SolveBallisticAngle(float v, float g, float x, float y,
            out float angle1, out float angle2)
        {
            // Solve: y = x*tan(θ) - (g*x²)/(2*v²*cos²(θ))
            // Using quadratic formula after substitution

            float v2 = v * v;
            float v4 = v2 * v2;
            float gx = g * x;
            float discriminant = v4 - g * (g * x * x + 2 * y * v2);

            if (discriminant < 0)
            {
                angle1 = angle2 = 0;
                return false;
            }

            float sqrtD = MathF.Sqrt(discriminant);
            angle1 = MathF.Atan((v2 - sqrtD) / gx);
            angle2 = MathF.Atan((v2 + sqrtD) / gx);

            return true;
        }

        #endregion
    }

    /// <summary>
    /// Real-time bullet simulation for games.
    /// </summary>
    public class BulletSimulator
    {
        private readonly List<ActiveBullet> _bullets = new();
        private readonly BallisticsSystem _ballistics = new();

        /// <summary>Environment conditions.</summary>
        public BallisticEnvironment Environment { get; set; } = new();

        /// <summary>Gets active bullets.</summary>
        public IReadOnlyList<ActiveBullet> ActiveBullets => _bullets;

        /// <summary>Event fired when a bullet hits something.</summary>
        public event Action<ActiveBullet, Vector3, Vector3>? OnBulletHit;

        /// <summary>Collision check function (position, direction, maxDist) -> (hit, point, normal).</summary>
        public Func<Vector3, Vector3, float, (bool hit, Vector3 point, Vector3 normal)>? CollisionCheck { get; set; }

        /// <summary>
        /// Spawns a new bullet.
        /// </summary>
        public ActiveBullet SpawnBullet(
            Projectile projectile,
            Vector3 position,
            Vector3 direction)
        {
            var bullet = new ActiveBullet
            {
                Projectile = projectile,
                Position = position,
                Velocity = Vector3.Normalize(direction) * projectile.MuzzleVelocity,
                SpinRate = projectile.SpinRate,
                TimeAlive = 0,
                Distance = 0,
                IsActive = true
            };

            _bullets.Add(bullet);
            return bullet;
        }

        /// <summary>
        /// Updates all active bullets.
        /// </summary>
        public void Update(float deltaTime)
        {
            float speedOfSound = Environment.SpeedOfSound;

            for (int i = _bullets.Count - 1; i >= 0; i--)
            {
                var bullet = _bullets[i];
                if (!bullet.IsActive)
                {
                    _bullets.RemoveAt(i);
                    continue;
                }

                // Store previous position for collision detection
                Vector3 prevPosition = bullet.Position;

                // Calculate forces
                float speed = bullet.Velocity.Length();
                if (speed < 1f)
                {
                    bullet.IsActive = false;
                    continue;
                }

                // Drag
                float mach = speed / speedOfSound;
                float dragCoeff = bullet.Projectile.DragCoefficient;
                if (mach > 0.8f && mach < 1.2f)
                    dragCoeff *= 1 + 0.5f * ((mach - 0.8f) / 0.4f);

                float dragMagnitude = 0.5f * Environment.AirDensity * speed * speed *
                                     dragCoeff * bullet.Projectile.CrossSectionalArea;
                Vector3 dragForce = -Vector3.Normalize(bullet.Velocity) * dragMagnitude;

                // Wind effect
                Vector3 relativeVel = bullet.Velocity - Environment.Wind;
                float windDrag = 0.5f * Environment.AirDensity * relativeVel.LengthSquared() *
                                dragCoeff * bullet.Projectile.CrossSectionalArea;
                Vector3 windForce = -Vector3.Normalize(relativeVel) * windDrag - dragForce;

                // Integrate
                Vector3 acceleration = Environment.Gravity +
                                       (dragForce + windForce) / bullet.Projectile.Mass;
                bullet.Velocity += acceleration * deltaTime;
                bullet.Position += bullet.Velocity * deltaTime;
                bullet.TimeAlive += deltaTime;
                bullet.Distance += speed * deltaTime;

                // Collision check
                if (CollisionCheck != null)
                {
                    Vector3 dir = bullet.Position - prevPosition;
                    float dist = dir.Length();
                    if (dist > 0.001f)
                    {
                        var (hit, point, normal) = CollisionCheck(prevPosition, dir / dist, dist);
                        if (hit)
                        {
                            bullet.IsActive = false;
                            OnBulletHit?.Invoke(bullet, point, normal);
                        }
                    }
                }

                // Ground check
                if (bullet.Position.Y < 0)
                {
                    bullet.IsActive = false;
                    OnBulletHit?.Invoke(bullet, new Vector3(bullet.Position.X, 0, bullet.Position.Z), Vector3.UnitY);
                }

                // Timeout
                if (bullet.TimeAlive > 10f)
                {
                    bullet.IsActive = false;
                }
            }
        }

        /// <summary>
        /// Clears all bullets.
        /// </summary>
        public void Clear()
        {
            _bullets.Clear();
        }
    }

    /// <summary>
    /// Active bullet in simulation.
    /// </summary>
    public class ActiveBullet
    {
        /// <summary>Projectile configuration.</summary>
        public Projectile Projectile { get; set; } = new();

        /// <summary>Current position.</summary>
        public Vector3 Position { get; set; }

        /// <summary>Current velocity.</summary>
        public Vector3 Velocity { get; set; }

        /// <summary>Current spin rate.</summary>
        public float SpinRate { get; set; }

        /// <summary>Time since spawn.</summary>
        public float TimeAlive { get; set; }

        /// <summary>Distance traveled.</summary>
        public float Distance { get; set; }

        /// <summary>Whether bullet is still active.</summary>
        public bool IsActive { get; set; } = true;

        /// <summary>Current speed.</summary>
        public float Speed => Velocity.Length();

        /// <summary>Current kinetic energy.</summary>
        public float KineticEnergy => 0.5f * Projectile.Mass * Velocity.LengthSquared();

        /// <summary>Custom user data.</summary>
        public object? UserData { get; set; }
    }

    /// <summary>
    /// Penetration calculator for terminal ballistics.
    /// </summary>
    public static class TerminalBallistics
    {
        /// <summary>
        /// Calculates penetration depth in a material.
        /// </summary>
        /// <param name="bullet">The projectile.</param>
        /// <param name="impactVelocity">Impact velocity in m/s.</param>
        /// <param name="targetHardness">Target hardness (BHN, 0-500).</param>
        /// <param name="targetDensity">Target density kg/m³.</param>
        /// <returns>Penetration depth in meters.</returns>
        public static float CalculatePenetration(
            Projectile bullet,
            float impactVelocity,
            float targetHardness = 150f,
            float targetDensity = 7800f)
        {
            // Simplified penetration model based on de Marre formula
            float energy = 0.5f * bullet.Mass * impactVelocity * impactVelocity;
            float area = bullet.CrossSectionalArea;

            // Penetration roughly proportional to energy/area and inversely to hardness
            float penetration = (energy / area) / (targetHardness * targetDensity * 0.001f);

            return MathF.Max(0, penetration);
        }

        /// <summary>
        /// Calculates bullet expansion ratio.
        /// </summary>
        /// <param name="impactVelocity">Impact velocity in m/s.</param>
        /// <param name="isHollowPoint">Whether bullet is hollow point.</param>
        /// <returns>Expansion ratio (1.0 = no expansion).</returns>
        public static float CalculateExpansion(float impactVelocity, bool isHollowPoint)
        {
            if (!isHollowPoint) return 1.0f;

            // Hollow points expand based on velocity
            if (impactVelocity < 250) return 1.0f;
            if (impactVelocity > 500) return 2.0f;

            return 1.0f + (impactVelocity - 250) / 250;
        }

        /// <summary>
        /// Calculates wound cavity volume (for simulation purposes).
        /// </summary>
        /// <param name="bullet">The projectile.</param>
        /// <param name="impactVelocity">Impact velocity.</param>
        /// <param name="penetration">Penetration depth.</param>
        /// <returns>Cavity volume in cm³.</returns>
        public static float CalculateCavityVolume(
            Projectile bullet,
            float impactVelocity,
            float penetration)
        {
            float energy = 0.5f * bullet.Mass * impactVelocity * impactVelocity;
            float diameter = bullet.Diameter * 100; // to cm

            // Simplified temporary cavity calculation
            return MathF.PI * diameter * diameter * penetration * 100 * (energy / 1000);
        }

        /// <summary>
        /// Estimates ricochet angle threshold.
        /// </summary>
        /// <param name="targetHardness">Surface hardness BHN.</param>
        /// <param name="bulletHardness">Bullet hardness BHN.</param>
        /// <returns>Maximum angle in degrees for ricochet.</returns>
        public static float CalculateRicochetAngle(float targetHardness, float bulletHardness = 150)
        {
            // Harder surfaces cause more ricochets at higher angles
            float ratio = targetHardness / bulletHardness;
            return 5 + 20 * MathF.Min(ratio, 3);
        }
    }
}
