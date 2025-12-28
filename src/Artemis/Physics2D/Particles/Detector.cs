using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Artemis.Physics2D.Particles
{
    /// <summary>
    /// Types of detectors in a particle collider (based on LHC).
    /// </summary>
    public enum DetectorType2D
    {
        Tracker,           // Tracks particle paths (silicon pixels)
        Calorimeter,       // Measures energy (electromagnetic + hadronic)
        MuonChamber,       // Detects muons (outer layer)
        InteractionPoint,  // Where collisions happen
        TriggerSystem,     // Fast event selection
        VertexDetector     // Precise vertex reconstruction
    }

    /// <summary>
    /// Detector that records particles passing through.
    /// Similar to ATLAS or CMS detectors at LHC.
    /// Thread-safe for parallel detection.
    /// </summary>
    public class Detector2D
    {
        public DetectorType2D Type { get; set; }
        public string Name { get; set; }
        public Vector2D Position { get; set; }
        public double Radius { get; set; }
        public uint Color { get; set; }

        // Trigger zone body (for physics world integration)
        public RigidBody2D? TriggerZone { get; private set; }

        // Detection data (thread-safe)
        public ConcurrentBag<ParticleDetection2D> Detections { get; } = new();
        private double _totalEnergyDetected;
        private int _particleCount;
        private readonly HashSet<int> _detectedParticleIds = new();
        private readonly object _detectionLock = new();

        public double TotalEnergyDetected => _totalEnergyDetected;
        public int ParticleCount => _particleCount;

        // Detector efficiency (0-1)
        public double Efficiency { get; set; } = 0.95;

        // Energy resolution (relative uncertainty)
        public double EnergyResolution { get; set; } = 0.03; // 3%

        // Position resolution (meters)
        public double PositionResolution { get; set; } = 0.001; // 1mm

        public Detector2D(DetectorType2D type, string name, Vector2D position, double radius)
        {
            Type = type;
            Name = name;
            Position = position;
            Radius = radius;

            Color = type switch
            {
                DetectorType2D.Tracker => 0x00FFFFFF,        // Cyan
                DetectorType2D.Calorimeter => 0xFFFF00FF,    // Yellow
                DetectorType2D.MuonChamber => 0xFF00FFFF,    // Magenta
                DetectorType2D.InteractionPoint => 0xFF0000FF, // Red
                DetectorType2D.TriggerSystem => 0x00FF00FF,  // Green
                DetectorType2D.VertexDetector => 0xFFA500FF, // Orange
                _ => 0xFFFFFFFF // White
            };

            SetupTriggerZone();
        }

        private void SetupTriggerZone()
        {
            // Create trigger zone (sensor, no physics response)
            TriggerZone = new RigidBody2D(Position)
            {
                BodyType = BodyType2D.Static,
                Shape = new CircleShape(Radius),
                IsActive = true,
                UserData = this
            };
            TriggerZone.CategoryBits = 0x0002; // Detector category
            TriggerZone.MaskBits = 0x0001;     // Detect particles
        }

        /// <summary>
        /// Check if particle is within detector volume.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool ContainsParticle(Particle2D particle)
        {
            double distance = (particle.Body.Position - Position).Magnitude;
            return distance <= Radius;
        }

        /// <summary>
        /// Record a particle detection (thread-safe).
        /// </summary>
        public void RecordDetection(Particle2D particle, double simulationTime)
        {
            // Check efficiency (stochastic)
            if (Random.Shared.NextDouble() > Efficiency)
                return;

            // Avoid double-counting the same particle
            lock (_detectionLock)
            {
                if (_detectedParticleIds.Contains(particle.Id))
                    return;
                _detectedParticleIds.Add(particle.Id);
            }

            // Apply measurement uncertainties
            double measuredEnergy = ApplyEnergyResolution(particle.GetKineticEnergy());
            Vector2D measuredPosition = ApplyPositionResolution(particle.Body.Position);
            Vector2D measuredVelocity = ApplyVelocityResolution(particle.Body.Velocity);

            var detection = new ParticleDetection2D
            {
                Time = simulationTime,
                ParticleType = particle.Properties.Type,
                ParticleId = particle.Id,
                Position = measuredPosition,
                Velocity = measuredVelocity,
                Energy = measuredEnergy,
                Momentum = particle.GetMomentum(),
                Charge = particle.Properties.Charge,
                DetectorName = Name,
                DetectorType = Type
            };

            Detections.Add(detection);
            Interlocked.Add(ref _particleCount, 1);

            // Thread-safe energy update
            double oldValue, newValue;
            do
            {
                oldValue = _totalEnergyDetected;
                newValue = oldValue + measuredEnergy;
            } while (Interlocked.CompareExchange(
                ref _totalEnergyDetected, newValue, oldValue) != oldValue);

            // Update particle's energy deposit tracking
            particle.TotalEnergyDeposited += measuredEnergy;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double ApplyEnergyResolution(double trueEnergy)
        {
            // Gaussian smearing
            double sigma = trueEnergy * EnergyResolution;
            return trueEnergy + sigma * RandomGaussian();
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private Vector2D ApplyPositionResolution(Vector2D truePosition)
        {
            return new Vector2D(
                truePosition.X + PositionResolution * RandomGaussian(),
                truePosition.Y + PositionResolution * RandomGaussian()
            );
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private Vector2D ApplyVelocityResolution(Vector2D trueVelocity)
        {
            // Velocity resolution derived from position resolution
            double sigma = trueVelocity.Magnitude * EnergyResolution * 0.5;
            return new Vector2D(
                trueVelocity.X + sigma * RandomGaussian(),
                trueVelocity.Y + sigma * RandomGaussian()
            );
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double RandomGaussian()
        {
            // Box-Muller transform
            double u1 = 1.0 - Random.Shared.NextDouble();
            double u2 = 1.0 - Random.Shared.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
        }

        /// <summary>
        /// Clear all detections.
        /// </summary>
        public void Clear()
        {
            while (Detections.TryTake(out _)) { }
            _totalEnergyDetected = 0;
            _particleCount = 0;
            lock (_detectionLock)
            {
                _detectedParticleIds.Clear();
            }
        }

        /// <summary>
        /// Reset per-event tracking (call between events).
        /// </summary>
        public void ResetEvent()
        {
            lock (_detectionLock)
            {
                _detectedParticleIds.Clear();
            }
        }

        public string GetSummary()
        {
            if (_particleCount == 0)
                return $"{Name}: No detections";

            return $"{Name}: {_particleCount} particles, {_totalEnergyDetected:F1} GeV";
        }

        /// <summary>
        /// Get detections of a specific particle type.
        /// </summary>
        public IEnumerable<ParticleDetection2D> GetDetectionsByType(ParticleType type)
        {
            foreach (var d in Detections)
            {
                if (d.ParticleType == type)
                    yield return d;
            }
        }

        /// <summary>
        /// Calculate invariant mass from two detected particles.
        /// </summary>
        public static double CalculateInvariantMass(ParticleDetection2D p1, ParticleDetection2D p2)
        {
            // E² = (E1 + E2)² - (p1 + p2)²c²
            // Simplified for c=1 units
            double e1 = p1.Energy;
            double e2 = p2.Energy;
            double px = p1.Momentum * Math.Cos(p1.Velocity.Angle) + p2.Momentum * Math.Cos(p2.Velocity.Angle);
            double py = p1.Momentum * Math.Sin(p1.Velocity.Angle) + p2.Momentum * Math.Sin(p2.Velocity.Angle);
            double p2Total = px * px + py * py;

            double m2 = (e1 + e2) * (e1 + e2) - p2Total;
            return m2 > 0 ? Math.Sqrt(m2) : 0;
        }
    }

    /// <summary>
    /// Record of a particle detection.
    /// </summary>
    public class ParticleDetection2D
    {
        public double Time { get; set; }
        public ParticleType ParticleType { get; set; }
        public int ParticleId { get; set; }
        public Vector2D Position { get; set; }
        public Vector2D Velocity { get; set; }
        public double Energy { get; set; }
        public double Momentum { get; set; }
        public int Charge { get; set; }
        public string DetectorName { get; set; } = "";
        public DetectorType2D DetectorType { get; set; }

        public override string ToString()
            => $"[{Time:F3}s] {ParticleType} @ {Position}, E={Energy:F2} GeV";
    }

    /// <summary>
    /// Collision event between particles.
    /// Represents a high-energy collision like those at the LHC.
    /// </summary>
    public class ParticleCollisionEvent2D
    {
        public double Time { get; set; }
        public Particle2D Particle1 { get; set; }
        public Particle2D Particle2 { get; set; }
        public Vector2D CollisionPoint { get; set; }
        public double TotalEnergyBefore { get; set; }
        public double TotalEnergyAfter { get; set; }
        public double CenterOfMassEnergy { get; set; }
        public List<Particle2D> Products { get; } = new();

        // Event classification
        public string EventType { get; set; } = "Unknown";
        public double InvariantMass { get; set; }

        public ParticleCollisionEvent2D(Particle2D p1, Particle2D p2, Vector2D point, double time)
        {
            Time = time;
            Particle1 = p1;
            Particle2 = p2;
            CollisionPoint = point;

            // Calculate center-of-mass energy (like √s at LHC)
            double e1 = p1.GetKineticEnergy();
            double e2 = p2.GetKineticEnergy();
            TotalEnergyBefore = e1 + e2;

            // Simplified COM energy calculation
            // √s = √(2 * m1 * m2 * c⁴ + 2 * E1 * E2 - 2 * p1·p2 * c²)
            var p1Vec = p1.GetMomentumVector();
            var p2Vec = p2.GetMomentumVector();
            double m1 = p1.Properties.Mass;
            double m2 = p2.Properties.Mass;

            CenterOfMassEnergy = Math.Sqrt(
                2 * m1 * m2 +
                2 * e1 * e2 -
                2 * Vector2D.Dot(p1Vec, p2Vec)
            );

            // Calculate invariant mass
            InvariantMass = CalculateInvariantMass(p1, p2);

            // Classify event
            ClassifyEvent();
        }

        private double CalculateInvariantMass(Particle2D p1, Particle2D p2)
        {
            var totalP = p1.GetMomentumVector() + p2.GetMomentumVector();
            double totalE = p1.GetKineticEnergy() + p2.GetKineticEnergy() +
                           p1.Properties.Mass + p2.Properties.Mass;

            double m2 = totalE * totalE - totalP.MagnitudeSquared;
            return m2 > 0 ? Math.Sqrt(m2) : 0;
        }

        private void ClassifyEvent()
        {
            // Simple classification based on particle types
            var t1 = Particle1.Properties.Type;
            var t2 = Particle2.Properties.Type;

            if (t1 == ParticleType.Proton && t2 == ParticleType.Proton)
            {
                EventType = "pp collision";
            }
            else if ((t1 == ParticleType.Electron && t2 == ParticleType.Positron) ||
                     (t1 == ParticleType.Positron && t2 == ParticleType.Electron))
            {
                EventType = "e+e- annihilation";
            }
            else if (t1 == ParticleType.Photon || t2 == ParticleType.Photon)
            {
                EventType = "Photon interaction";
            }
            else
            {
                EventType = $"{t1}-{t2} collision";
            }
        }

        public void AnalyzeProducts(IEnumerable<Particle2D> allParticles)
        {
            Products.Clear();

            // Find particles created near collision point recently
            foreach (var p in allParticles)
            {
                if ((p.Body.Position - CollisionPoint).Magnitude < 2.0 &&
                    p.Age < 0.1 &&
                    p != Particle1 && p != Particle2)
                {
                    Products.Add(p);
                }
            }

            TotalEnergyAfter = 0;
            foreach (var p in Products)
            {
                TotalEnergyAfter += p.GetKineticEnergy();
            }
        }

        public string GetSummary()
        {
            string p1Type = Particle1.Properties.Symbol;
            string p2Type = Particle2.Properties.Symbol;

            return $"{p1Type} + {p2Type} → {Products.Count} products, " +
                   $"√s = {CenterOfMassEnergy:F1} GeV, " +
                   $"M_inv = {InvariantMass:F1} GeV, " +
                   $"E_in = {TotalEnergyBefore:F1} GeV";
        }
    }

    /// <summary>
    /// Energy tracking and conservation for scientific simulations.
    /// </summary>
    public class EnergyTracker2D
    {
        public double KineticEnergy { get; private set; }
        public double PotentialEnergy { get; private set; }
        public double TotalEnergy => KineticEnergy + PotentialEnergy;
        public double RotationalEnergy { get; private set; }

        public Vector2D GravityDirection { get; set; } = new Vector2D(0, -1);
        public double GravityMagnitude { get; set; } = 9.81;

        private double _initialEnergy;
        private bool _initialized;

        public void UpdateEnergy(IEnumerable<RigidBody2D> bodies)
        {
            KineticEnergy = 0;
            PotentialEnergy = 0;
            RotationalEnergy = 0;

            foreach (var body in bodies)
            {
                if (body.BodyType == BodyType2D.Static)
                    continue;

                // Kinetic energy: KE = 0.5 * m * v²
                KineticEnergy += 0.5 * body.Mass * body.Velocity.MagnitudeSquared;

                // Rotational kinetic energy: RE = 0.5 * I * ω²
                RotationalEnergy += 0.5 * body.Inertia * body.AngularVelocity * body.AngularVelocity;

                // Potential energy: PE = m * g * h
                double height = Vector2D.Dot(body.Position, -GravityDirection);
                PotentialEnergy += body.Mass * GravityMagnitude * height;
            }

            if (!_initialized)
            {
                _initialEnergy = TotalEnergy;
                _initialized = true;
            }
        }

        public double GetEnergyDrift()
        {
            return TotalEnergy - _initialEnergy;
        }

        public double GetEnergyDriftPercentage()
        {
            if (_initialEnergy == 0)
                return 0;

            return (GetEnergyDrift() / _initialEnergy) * 100;
        }

        public void Reset()
        {
            _initialized = false;
            _initialEnergy = 0;
        }
    }

    /// <summary>
    /// Data export for scientific analysis.
    /// </summary>
    public class SimulationDataExporter2D
    {
        public class FrameData
        {
            public double Time { get; set; }
            public List<ParticleState> Particles { get; set; } = new();
        }

        public class ParticleState
        {
            public int Id { get; set; }
            public ParticleType Type { get; set; }
            public Vector2D Position { get; set; }
            public Vector2D Velocity { get; set; }
            public double Energy { get; set; }
            public double Momentum { get; set; }
            public int Charge { get; set; }
        }

        private List<FrameData> _frames = new();
        private double _currentTime;

        public void RecordFrame(IEnumerable<Particle2D> particles, double deltaTime)
        {
            _currentTime += deltaTime;

            var frame = new FrameData { Time = _currentTime };

            foreach (var p in particles)
            {
                frame.Particles.Add(new ParticleState
                {
                    Id = p.Id,
                    Type = p.Properties.Type,
                    Position = p.Body.Position,
                    Velocity = p.Body.Velocity,
                    Energy = p.GetKineticEnergy(),
                    Momentum = p.GetMomentum(),
                    Charge = p.Properties.Charge
                });
            }

            _frames.Add(frame);
        }

        public string ExportToCSV()
        {
            var csv = new System.Text.StringBuilder();
            csv.AppendLine("Time,ParticleId,Type,PosX,PosY,VelX,VelY,Energy,Momentum,Charge");

            foreach (var frame in _frames)
            {
                foreach (var p in frame.Particles)
                {
                    csv.AppendLine($"{frame.Time:F6},{p.Id},{p.Type},{p.Position.X:F4},{p.Position.Y:F4}," +
                                  $"{p.Velocity.X:F4},{p.Velocity.Y:F4},{p.Energy:F4},{p.Momentum:F4},{p.Charge}");
                }
            }

            return csv.ToString();
        }

        public void Clear()
        {
            _frames.Clear();
            _currentTime = 0;
        }

        public int FrameCount => _frames.Count;
    }
}
