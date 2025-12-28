using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Artemis.Physics2D.Particles
{
    /// <summary>
    /// The particle accelerator ring (like the LHC's 27km ring).
    /// Uses magnetic field effectors to keep particles in orbit.
    /// Multi-threaded for high particle counts.
    /// </summary>
    public class Accelerator2D
    {
        public Vector2D Center { get; set; }
        public double Radius { get; set; }
        public double BeamEnergy { get; set; } // Energy in GeV

        // Two counter-rotating beams (like LHC)
        public ParticleBeam Beam1 { get; } = new() { Name = "Beam 1", IsClockwise = true };
        public ParticleBeam Beam2 { get; } = new() { Name = "Beam 2", IsClockwise = false };

        // Magnetic confinement
        private RadialForceEffector2D _magneticField;
        private VortexEffector2D _clockwiseVortex;
        private VortexEffector2D _counterClockwiseVortex;

        // RF cavities (acceleration)
        private List<RFCavity2D> _rfCavities = new();

        // Focusing quadrupoles
        private double _focusingStrength = 10.0;

        public Accelerator2D(Vector2D center, double radius, double beamEnergy)
        {
            Center = center;
            Radius = radius;
            BeamEnergy = beamEnergy;

            // Create magnetic field to keep particles in orbit (centripetal)
            _magneticField = new RadialForceEffector2D(
                center: center,
                radius: radius + 5,
                strength: -50.0, // Attractive force toward center
                falloff: RadialFalloff2D.Linear
            );

            // Vortex for beam 1 (clockwise)
            _clockwiseVortex = new VortexEffector2D(
                center: center,
                radius: radius + 2,
                strength: 20.0,
                clockwise: true
            );

            // Vortex for beam 2 (counter-clockwise)
            _counterClockwiseVortex = new VortexEffector2D(
                center: center,
                radius: radius + 2,
                strength: 20.0,
                clockwise: false
            );

            // Create RF cavities around the ring
            SetupRFCavities(8);
        }

        private void SetupRFCavities(int count)
        {
            _rfCavities.Clear();
            for (int i = 0; i < count; i++)
            {
                double angle = i * (2 * Math.PI / count);
                var position = Center + Vector2D.FromAngle(angle, Radius);
                _rfCavities.Add(new RFCavity2D(position, 2.0, BeamEnergy * 0.01));
            }
        }

        public RadialForceEffector2D GetMagneticField() => _magneticField;
        public VortexEffector2D GetClockwiseVortex() => _clockwiseVortex;
        public VortexEffector2D GetCounterClockwiseVortex() => _counterClockwiseVortex;

        /// <summary>
        /// Register effectors with a physics world.
        /// </summary>
        public void RegisterEffectors(PhysicsWorld2D world)
        {
            world.AddAreaEffector(_magneticField);
            // Note: Vortex effectors are applied selectively based on beam direction
        }

        /// <summary>
        /// Inject particles into beam 1 (clockwise).
        /// </summary>
        public void InjectBeam1(ParticleType type, int count)
        {
            for (int i = 0; i < count; i++)
            {
                double angle = (double)i / count * Math.PI * 2;
                Vector2D position = Center + Vector2D.FromAngle(angle, Radius);

                // Velocity tangent to circle (clockwise)
                Vector2D tangent = Vector2D.FromAngle(angle - Math.PI / 2);
                double speed = CalculateBeamSpeed(type);
                Vector2D velocity = tangent * speed;

                var particle = new Particle2D(type, position, velocity);
                Beam1.Add(particle);
            }
        }

        /// <summary>
        /// Inject particles into beam 2 (counter-clockwise).
        /// </summary>
        public void InjectBeam2(ParticleType type, int count)
        {
            for (int i = 0; i < count; i++)
            {
                double angle = (double)i / count * Math.PI * 2;
                Vector2D position = Center + Vector2D.FromAngle(angle, Radius);

                // Velocity tangent to circle (counter-clockwise)
                Vector2D tangent = Vector2D.FromAngle(angle + Math.PI / 2);
                double speed = CalculateBeamSpeed(type);
                Vector2D velocity = tangent * speed;

                var particle = new Particle2D(type, position, velocity);
                Beam2.Add(particle);
            }
        }

        /// <summary>
        /// Calculate particle speed based on beam energy.
        /// E = 0.5 * m * v^2, so v = sqrt(2*E/m)
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double CalculateBeamSpeed(ParticleType type)
        {
            var props = ParticleProperties.Create(type);
            // Avoid division by zero for massless particles
            double mass = Math.Max(props.Mass, 0.0001);
            return Math.Sqrt(2 * BeamEnergy / mass);
        }

        /// <summary>
        /// Get interaction points where beams should cross.
        /// At LHC, there are 4 main interaction points (ATLAS, CMS, ALICE, LHCb).
        /// </summary>
        public List<Vector2D> GetInteractionPoints()
        {
            return new List<Vector2D>
            {
                Center + new Vector2D(Radius, 0),      // IP1 (ATLAS-like)
                Center + new Vector2D(0, Radius),      // IP2 (ALICE-like)
                Center + new Vector2D(-Radius, 0),     // IP3 (CMS-like)
                Center + new Vector2D(0, -Radius)      // IP4 (LHCb-like)
            };
        }

        /// <summary>
        /// Apply focusing to keep particles in beam (parallel, simplified quadrupole model).
        /// </summary>
        public void ApplyBeamFocusingParallel(double deltaTime)
        {
            var allParticles = GetAllParticles().ToArray();

            Parallel.ForEach(allParticles, particle =>
            {
                ApplyFocusing(particle, deltaTime);
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void ApplyFocusing(Particle2D particle, double deltaTime)
        {
            // Calculate distance from ideal orbit
            double distanceFromCenter = (particle.Body.Position - Center).Magnitude;
            double deviation = distanceFromCenter - Radius;

            if (Math.Abs(deviation) > 0.1)
            {
                // Apply restoring force toward ideal orbit
                Vector2D toCenter = (Center - particle.Body.Position).Normalized;
                Vector2D restoringForce = toCenter * (deviation * _focusingStrength);
                particle.Body.ApplyForce(restoringForce * deltaTime);
            }
        }

        /// <summary>
        /// Apply RF cavity acceleration (parallel).
        /// </summary>
        public void ApplyRFAccelerationParallel(double deltaTime, double phase = 0)
        {
            var allParticles = GetAllParticles().ToArray();

            Parallel.ForEach(allParticles, particle =>
            {
                foreach (var cavity in _rfCavities)
                {
                    cavity.ApplyAcceleration(particle, deltaTime, phase);
                }
            });
        }

        /// <summary>
        /// Update all particles (parallel).
        /// </summary>
        public void UpdateParticlesParallel(double deltaTime)
        {
            // Update Beam 1
            Parallel.ForEach(Beam1.Particles, p => p.Update(deltaTime));

            // Update Beam 2
            Parallel.ForEach(Beam2.Particles, p => p.Update(deltaTime));

            // Remove decayed particles
            Beam1.RemoveDecayed();
            Beam2.RemoveDecayed();
        }

        /// <summary>
        /// Get all particles in the accelerator.
        /// </summary>
        public IEnumerable<Particle2D> GetAllParticles()
        {
            return Beam1.Particles.Concat(Beam2.Particles);
        }

        /// <summary>
        /// Clear all beams.
        /// </summary>
        public void ClearBeams()
        {
            Beam1.Clear();
            Beam2.Clear();
        }

        /// <summary>
        /// Get beam statistics.
        /// </summary>
        public (int beam1Count, int beam2Count, double avgEnergy, double luminosity) GetStatistics()
        {
            int b1 = Beam1.Count;
            int b2 = Beam2.Count;
            double avgE = 0;

            var all = GetAllParticles().ToList();
            if (all.Count > 0)
            {
                avgE = all.Average(p => p.GetKineticEnergy());
            }

            // Simplified luminosity calculation
            double luminosity = (b1 * b2) / (Math.PI * Radius * Radius) * 1e30;

            return (b1, b2, avgE, luminosity);
        }

        /// <summary>
        /// Set focusing strength.
        /// </summary>
        public void SetFocusingStrength(double strength)
        {
            _focusingStrength = strength;
        }

        /// <summary>
        /// Set magnetic field strength.
        /// </summary>
        public void SetMagneticFieldStrength(double strength)
        {
            _magneticField.Strength = strength;
        }
    }

    /// <summary>
    /// RF (Radio Frequency) cavity for particle acceleration.
    /// </summary>
    public class RFCavity2D
    {
        public Vector2D Position { get; set; }
        public double Radius { get; set; }
        public double AccelerationStrength { get; set; }
        public double Frequency { get; set; } = 400e6; // 400 MHz like LHC

        public RFCavity2D(Vector2D position, double radius, double strength)
        {
            Position = position;
            Radius = radius;
            AccelerationStrength = strength;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void ApplyAcceleration(Particle2D particle, double deltaTime, double phase = 0)
        {
            var delta = particle.Body.Position - Position;
            double distance = delta.Magnitude;

            if (distance > Radius)
                return;

            // Only accelerate charged particles
            if (particle.Properties.Charge == 0)
                return;

            // Apply acceleration in the direction of motion
            var velocity = particle.Body.Velocity;
            if (velocity.MagnitudeSquared < 0.001)
                return;

            var direction = velocity.Normalized;
            double strength = AccelerationStrength * particle.Properties.Charge;

            // RF phase modulation
            double phaseModulation = Math.Sin(phase);
            particle.Body.ApplyForce(direction * strength * phaseModulation * deltaTime);
        }
    }

    /// <summary>
    /// Experiment coordinator - manages the collision experiment.
    /// Multi-threaded collision detection and event recording.
    /// </summary>
    public class CollisionExperiment2D
    {
        public string Name { get; set; }
        public Accelerator2D Accelerator { get; }
        public List<Detector2D> Detectors { get; } = new();
        public ConcurrentBag<ParticleCollisionEvent2D> CollisionEvents { get; } = new();

        public double SimulationTime { get; private set; }
        public int TotalCollisions { get; private set; }

        // Data tracking
        private EnergyTracker2D _energyTracker = new();

        public CollisionExperiment2D(string name, Accelerator2D accelerator)
        {
            Name = name;
            Accelerator = accelerator;
            SetupDetectors();
        }

        private void SetupDetectors()
        {
            var interactionPoints = Accelerator.GetInteractionPoints();
            string[] detectorNames = { "ATLAS", "ALICE", "CMS", "LHCb" };

            for (int i = 0; i < interactionPoints.Count; i++)
            {
                var ip = interactionPoints[i];
                string name = i < detectorNames.Length ? detectorNames[i] : $"IP{i + 1}";

                // Inner tracker
                Detectors.Add(new Detector2D(
                    DetectorType2D.Tracker,
                    $"{name}_Tracker",
                    ip,
                    3.0
                ));

                // Calorimeter (outer layer)
                Detectors.Add(new Detector2D(
                    DetectorType2D.Calorimeter,
                    $"{name}_Calorimeter",
                    ip,
                    5.0
                ));

                // Interaction point marker
                Detectors.Add(new Detector2D(
                    DetectorType2D.InteractionPoint,
                    name,
                    ip,
                    0.5
                ));
            }
        }

        /// <summary>
        /// Update experiment (parallel).
        /// </summary>
        public void Update(double deltaTime, PhysicsWorld2D world)
        {
            SimulationTime += deltaTime;

            // Update energy tracking
            var allBodies = world.Bodies
                .Where(b => b.UserData is Particle2D)
                .ToList();
            _energyTracker.UpdateEnergy(allBodies);

            // Detect particles in detectors (parallel)
            DetectParticlesParallel();
        }

        private void DetectParticlesParallel()
        {
            var allParticles = Accelerator.GetAllParticles().ToArray();

            Parallel.ForEach(Detectors, detector =>
            {
                foreach (var particle in allParticles)
                {
                    if (detector.ContainsParticle(particle))
                    {
                        detector.RecordDetection(particle, SimulationTime);
                    }
                }
            });
        }

        public void RecordCollision(Particle2D p1, Particle2D p2, Vector2D point)
        {
            var collisionEvent = new ParticleCollisionEvent2D(p1, p2, point, SimulationTime);
            CollisionEvents.Add(collisionEvent);
            Interlocked.Increment(ref _totalCollisions);
        }

        private int _totalCollisions;

        public string GetEnergyReport()
        {
            return $"KE: {_energyTracker.KineticEnergy:F1} GeV | " +
                   $"Total: {_energyTracker.TotalEnergy:F1} GeV";
        }

        public string GetStatisticsReport()
        {
            var stats = Accelerator.GetStatistics();
            return $"Beam1: {stats.beam1Count} | Beam2: {stats.beam2Count} | " +
                   $"Avg E: {stats.avgEnergy:F1} GeV | L: {stats.luminosity:E2} cm⁻²s⁻¹";
        }
    }
}
