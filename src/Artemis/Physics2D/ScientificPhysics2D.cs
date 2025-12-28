using System;
using System.Collections.Generic;
using System.Text;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Advanced integration methods for scientific simulations.
    /// More accurate and stable than basic Euler integration.
    /// </summary>
    public static class AdvancedIntegration2D
    {
        /// <summary>
        /// Verlet integration - better energy conservation than Euler.
        /// Commonly used in molecular dynamics and particle systems.
        /// </summary>
        public static void VerletIntegration(RigidBody2D body, Vector2D acceleration, double deltaTime)
        {
            if (body.BodyType != BodyType2D.Dynamic || body.IsSleeping)
                return;

            // Store previous position if not available (first frame)
            Vector2D previousPosition = body.Position - body.Velocity * deltaTime;

            // Verlet integration: x(t+dt) = 2*x(t) - x(t-dt) + a*dt^2
            Vector2D newPosition = 2 * body.Position - previousPosition + acceleration * deltaTime * deltaTime;

            // Calculate velocity from position difference
            body.Velocity = (newPosition - body.Position) / deltaTime;
            body.Position = newPosition;
        }

        /// <summary>
        /// Runge-Kutta 4th order (RK4) - highly accurate integration.
        /// Used in aerospace, orbital mechanics, and precision simulations.
        /// </summary>
        public static void RK4Integration(RigidBody2D body, Func<Vector2D, Vector2D, Vector2D> accelerationFunc, double deltaTime)
        {
            if (body.BodyType != BodyType2D.Dynamic || body.IsSleeping)
                return;

            Vector2D x = body.Position;
            Vector2D v = body.Velocity;

            // RK4 steps
            Vector2D k1v = accelerationFunc(x, v);
            Vector2D k1x = v;

            Vector2D k2v = accelerationFunc(x + k1x * (deltaTime / 2), v + k1v * (deltaTime / 2));
            Vector2D k2x = v + k1v * (deltaTime / 2);

            Vector2D k3v = accelerationFunc(x + k2x * (deltaTime / 2), v + k2v * (deltaTime / 2));
            Vector2D k3x = v + k2v * (deltaTime / 2);

            Vector2D k4v = accelerationFunc(x + k3x * deltaTime, v + k3v * deltaTime);
            Vector2D k4x = v + k3v * deltaTime;

            // Weighted average
            body.Velocity += (k1v + 2 * k2v + 2 * k3v + k4v) * (deltaTime / 6);
            body.Position += (k1x + 2 * k2x + 2 * k3x + k4x) * (deltaTime / 6);
        }

        /// <summary>
        /// Semi-implicit Euler (Symplectic Euler) - better for game physics.
        /// More stable than explicit Euler, energy preserving.
        /// </summary>
        public static void SemiImplicitEuler(RigidBody2D body, Vector2D acceleration, double deltaTime)
        {
            if (body.BodyType != BodyType2D.Dynamic || body.IsSleeping)
                return;

            // Update velocity first
            body.Velocity += acceleration * deltaTime;

            // Then update position with new velocity
            body.Position += body.Velocity * deltaTime;
        }
    }

    /// <summary>
    /// Energy tracking and conservation for scientific simulations.
    /// Monitors kinetic, potential, and total energy of the system.
    /// </summary>
    public class EnergyTracker2D
    {
        public double KineticEnergy { get; private set; }
        public double PotentialEnergy { get; private set; }
        public double TotalEnergy => KineticEnergy + PotentialEnergy;
        public double RotationalEnergy { get; private set; }

        public Vector2D GravityDirection { get; set; } = new Vector2D(0, -1);
        public double GravityMagnitude { get; set; } = 9.81;

        public void UpdateEnergy(IEnumerable<RigidBody2D> bodies)
        {
            KineticEnergy = 0;
            PotentialEnergy = 0;
            RotationalEnergy = 0;

            foreach (var body in bodies)
            {
                if (body.BodyType != BodyType2D.Dynamic)
                    continue;

                // Kinetic energy: KE = 0.5 * m * v^2
                KineticEnergy += 0.5 * body.Mass * body.Velocity.MagnitudeSquared;

                // Rotational kinetic energy: RE = 0.5 * I * ω^2
                RotationalEnergy += 0.5 * body.Inertia * body.AngularVelocity * body.AngularVelocity;

                // Potential energy: PE = m * g * h
                double height = Vector2D.Dot(body.Position, -GravityDirection);
                PotentialEnergy += body.Mass * GravityMagnitude * height;
            }
        }

        public double GetEnergyDrift(double initialEnergy)
        {
            return TotalEnergy - initialEnergy;
        }

        public double GetEnergyDriftPercentage(double initialEnergy)
        {
            if (initialEnergy == 0)
                return 0;

            return (GetEnergyDrift(initialEnergy) / initialEnergy) * 100;
        }
    }

    /// <summary>
    /// Deterministic physics for reproducible simulations.
    /// Ensures same input always produces same output.
    /// </summary>
    public class DeterministicPhysics2D
    {
        private uint _seed;
        private Random _random;

        public DeterministicPhysics2D(uint seed = 12345)
        {
            _seed = seed;
            _random = new Random((int)seed);
        }

        public void Reset()
        {
            _random = new Random((int)_seed);
        }

        public double NextDouble()
        {
            return _random.NextDouble();
        }

        public double NextDouble(double min, double max)
        {
            return min + (max - min) * NextDouble();
        }

        public Vector2D NextVector2D(double minMag, double maxMag)
        {
            double angle = NextDouble() * 2 * Math.PI;
            double mag = NextDouble(minMag, maxMag);
            return new Vector2D(Math.Cos(angle) * mag, Math.Sin(angle) * mag);
        }
    }

    /// <summary>
    /// Data export for scientific analysis.
    /// Export simulation data to CSV or other formats.
    /// </summary>
    public class SimulationDataExporter2D
    {
        public class FrameData
        {
            public double Time { get; set; }
            public List<BodyState> Bodies { get; set; } = new();
        }

        public class BodyState
        {
            public string? Id { get; set; }
            public Vector2D Position { get; set; }
            public Vector2D Velocity { get; set; }
            public double Rotation { get; set; }
            public double AngularVelocity { get; set; }
            public double KineticEnergy { get; set; }
        }

        private List<FrameData> _frames = new();
        private double _currentTime = 0;

        public void RecordFrame(IEnumerable<RigidBody2D> bodies, double deltaTime)
        {
            _currentTime += deltaTime;

            var frame = new FrameData { Time = _currentTime };

            foreach (var body in bodies)
            {
                frame.Bodies.Add(new BodyState
                {
                    Id = body.Id,
                    Position = body.Position,
                    Velocity = body.Velocity,
                    Rotation = body.Rotation,
                    AngularVelocity = body.AngularVelocity,
                    KineticEnergy = 0.5 * body.Mass * body.Velocity.MagnitudeSquared +
                                   0.5 * body.Inertia * body.AngularVelocity * body.AngularVelocity
                });
            }

            _frames.Add(frame);
        }

        public string ExportToCSV()
        {
            var csv = new StringBuilder();
            csv.AppendLine("Time,BodyId,PosX,PosY,VelX,VelY,Rotation,AngularVel,KineticEnergy");

            foreach (var frame in _frames)
            {
                foreach (var body in frame.Bodies)
                {
                    csv.AppendLine($"{frame.Time},{body.Id},{body.Position.X},{body.Position.Y}," +
                                  $"{body.Velocity.X},{body.Velocity.Y},{body.Rotation}," +
                                  $"{body.AngularVelocity},{body.KineticEnergy}");
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

    /// <summary>
    /// Precision settings for scientific simulations.
    /// </summary>
    public class PrecisionSettings2D
    {
        public enum IntegrationMethod
        {
            Euler,
            SemiImplicitEuler,
            Verlet,
            RK4
        }

        public IntegrationMethod Method { get; set; } = IntegrationMethod.SemiImplicitEuler;
        public int ConstraintIterations { get; set; } = 10;
        public int VelocityIterations { get; set; } = 8;
        public int PositionIterations { get; set; } = 3;
        public double Tolerance { get; set; } = 0.0001;
    }
}
