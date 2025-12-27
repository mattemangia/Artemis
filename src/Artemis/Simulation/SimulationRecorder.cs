using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Artemis.Core;
using Artemis.Bodies;
using Artemis.Particles;

namespace Artemis.Simulation
{
    /// <summary>
    /// Records simulation state over time for analysis and playback.
    /// </summary>
    public class SimulationRecorder
    {
        #region Fields

        private readonly List<SimulationFrame> _frames;
        private readonly Dictionary<string, List<Vector3D>> _trajectories;
        private readonly Dictionary<string, BodyStatistics> _bodyStats;
        private double _currentTime;
        private bool _isRecording;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the recorded frames.
        /// </summary>
        public IReadOnlyList<SimulationFrame> Frames => _frames;

        /// <summary>
        /// Gets the total number of recorded frames.
        /// </summary>
        public int FrameCount => _frames.Count;

        /// <summary>
        /// Gets the total recorded time.
        /// </summary>
        public double TotalTime => _currentTime;

        /// <summary>
        /// Gets whether recording is active.
        /// </summary>
        public bool IsRecording => _isRecording;

        /// <summary>
        /// Gets or sets the recording interval (record every N updates).
        /// </summary>
        public int RecordInterval { get; set; } = 1;

        /// <summary>
        /// Gets or sets whether to record particle data.
        /// </summary>
        public bool RecordParticles { get; set; } = true;

        /// <summary>
        /// Gets or sets whether to track trajectories.
        /// </summary>
        public bool TrackTrajectories { get; set; } = true;

        /// <summary>
        /// Gets or sets the maximum frames to store (0 = unlimited).
        /// </summary>
        public int MaxFrames { get; set; } = 0;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new simulation recorder.
        /// </summary>
        public SimulationRecorder()
        {
            _frames = new List<SimulationFrame>();
            _trajectories = new Dictionary<string, List<Vector3D>>();
            _bodyStats = new Dictionary<string, BodyStatistics>();
            _currentTime = 0;
            _isRecording = false;
        }

        #endregion

        #region Recording

        /// <summary>
        /// Starts recording.
        /// </summary>
        public void StartRecording()
        {
            _isRecording = true;
        }

        /// <summary>
        /// Stops recording.
        /// </summary>
        public void StopRecording()
        {
            _isRecording = false;
        }

        /// <summary>
        /// Clears all recorded data.
        /// </summary>
        public void Clear()
        {
            _frames.Clear();
            _trajectories.Clear();
            _bodyStats.Clear();
            _currentTime = 0;
        }

        private int _updateCounter = 0;

        /// <summary>
        /// Records the current state of a physics world.
        /// </summary>
        public void RecordFrame(PhysicsWorld world, double deltaTime)
        {
            if (!_isRecording)
                return;

            _currentTime += deltaTime;
            _updateCounter++;

            if (_updateCounter % RecordInterval != 0)
                return;

            var frame = new SimulationFrame
            {
                Time = _currentTime,
                DeltaTime = deltaTime,
                FrameIndex = _frames.Count,
                BodyStates = new List<BodyState>()
            };

            // Record body states
            foreach (var body in world.Bodies)
            {
                var state = new BodyState
                {
                    Id = body.Id,
                    Position = body.Position,
                    Velocity = body.Velocity,
                    Rotation = body.Rotation,
                    AngularVelocity = body.AngularVelocity,
                    KineticEnergy = CalculateKineticEnergy(body),
                    PotentialEnergy = CalculatePotentialEnergy(body, world.Gravity)
                };
                frame.BodyStates.Add(state);

                // Track trajectory
                if (TrackTrajectories)
                {
                    if (!_trajectories.ContainsKey(body.Id))
                        _trajectories[body.Id] = new List<Vector3D>();
                    _trajectories[body.Id].Add(body.Position);
                }

                // Update statistics
                UpdateBodyStatistics(body);
            }

            // Calculate aggregate statistics
            frame.TotalKineticEnergy = frame.BodyStates.Sum(s => s.KineticEnergy);
            frame.TotalPotentialEnergy = frame.BodyStates.Sum(s => s.PotentialEnergy);
            frame.TotalEnergy = frame.TotalKineticEnergy + frame.TotalPotentialEnergy;

            // Enforce max frames limit
            if (MaxFrames > 0 && _frames.Count >= MaxFrames)
            {
                _frames.RemoveAt(0);
            }

            _frames.Add(frame);
        }

        /// <summary>
        /// Records the current state of a particle system.
        /// </summary>
        public void RecordFrame(ParticleSystem particles, double deltaTime)
        {
            if (!_isRecording || !RecordParticles)
                return;

            _currentTime += deltaTime;
            _updateCounter++;

            if (_updateCounter % RecordInterval != 0)
                return;

            var frame = new SimulationFrame
            {
                Time = _currentTime,
                DeltaTime = deltaTime,
                FrameIndex = _frames.Count,
                ParticleStates = new List<ParticleState>()
            };

            // Record particle states
            for (int i = 0; i < particles.MaxParticles; i++)
            {
                var p = particles.Particles[i];
                if (!p.IsAlive)
                    continue;

                frame.ParticleStates.Add(new ParticleState
                {
                    Index = i,
                    Position = p.Position,
                    Velocity = p.Velocity,
                    Mass = p.Mass,
                    Radius = p.Radius,
                    Lifetime = p.Lifetime,
                    Color = p.Color
                });
            }

            frame.ParticleCount = frame.ParticleStates.Count;

            if (MaxFrames > 0 && _frames.Count >= MaxFrames)
            {
                _frames.RemoveAt(0);
            }

            _frames.Add(frame);
        }

        private double CalculateKineticEnergy(IPhysicsBody body)
        {
            // KE = 0.5 * m * v²
            double linearKE = 0.5 * body.Mass * body.Velocity.MagnitudeSquared;

            // Rotational KE = 0.5 * I * ω² (simplified as scalar)
            double angularSpeed = body.AngularVelocity.Magnitude;
            double rotationalKE = 0.5 * body.Mass * angularSpeed * angularSpeed;

            return linearKE + rotationalKE;
        }

        private double CalculatePotentialEnergy(IPhysicsBody body, Vector3D gravity)
        {
            // PE = m * g * h (height relative to Y=0)
            double g = gravity.Magnitude;
            return body.Mass * g * body.Position.Y;
        }

        private void UpdateBodyStatistics(IPhysicsBody body)
        {
            if (!_bodyStats.ContainsKey(body.Id))
            {
                _bodyStats[body.Id] = new BodyStatistics { Id = body.Id };
            }

            var stats = _bodyStats[body.Id];
            double speed = body.Velocity.Magnitude;

            stats.SampleCount++;
            stats.MaxSpeed = Math.Max(stats.MaxSpeed, speed);
            stats.MinPosition = Vector3D.Min(stats.MinPosition, body.Position);
            stats.MaxPosition = Vector3D.Max(stats.MaxPosition, body.Position);

            // Running average
            stats.AverageSpeed = stats.AverageSpeed +
                (speed - stats.AverageSpeed) / stats.SampleCount;
        }

        #endregion

        #region Playback

        /// <summary>
        /// Gets a frame at a specific index.
        /// </summary>
        public SimulationFrame GetFrame(int index)
        {
            if (index < 0 || index >= _frames.Count)
                throw new ArgumentOutOfRangeException(nameof(index));
            return _frames[index];
        }

        /// <summary>
        /// Gets a frame at a specific time (interpolated).
        /// </summary>
        public SimulationFrame GetFrameAtTime(double time)
        {
            if (_frames.Count == 0)
                throw new InvalidOperationException("No frames recorded");

            // Find surrounding frames
            int index = _frames.FindIndex(f => f.Time >= time);

            if (index < 0)
                return _frames[_frames.Count - 1];
            if (index == 0)
                return _frames[0];

            // Could interpolate here if needed
            return _frames[index];
        }

        /// <summary>
        /// Gets the trajectory of a body.
        /// </summary>
        public IReadOnlyList<Vector3D> GetTrajectory(string bodyId)
        {
            return _trajectories.TryGetValue(bodyId, out var trajectory)
                ? trajectory
                : Array.Empty<Vector3D>();
        }

        /// <summary>
        /// Gets statistics for a body.
        /// </summary>
        public BodyStatistics? GetStatistics(string bodyId)
        {
            return _bodyStats.TryGetValue(bodyId, out var stats) ? stats : null;
        }

        /// <summary>
        /// Gets all body statistics.
        /// </summary>
        public IReadOnlyDictionary<string, BodyStatistics> GetAllStatistics()
        {
            return _bodyStats;
        }

        #endregion

        #region Analysis

        /// <summary>
        /// Calculates the total energy over time.
        /// </summary>
        public List<(double time, double energy)> GetEnergyOverTime()
        {
            return _frames
                .Select(f => (f.Time, f.TotalEnergy))
                .ToList();
        }

        /// <summary>
        /// Calculates energy conservation error (should be near zero for closed systems).
        /// </summary>
        public double CalculateEnergyConservationError()
        {
            if (_frames.Count < 2)
                return 0;

            double initialEnergy = _frames[0].TotalEnergy;
            double finalEnergy = _frames[_frames.Count - 1].TotalEnergy;

            if (Math.Abs(initialEnergy) < PhysicsConstants.Epsilon)
                return 0;

            return Math.Abs(finalEnergy - initialEnergy) / Math.Abs(initialEnergy);
        }

        /// <summary>
        /// Gets position over time for a specific body.
        /// </summary>
        public List<(double time, Vector3D position)> GetPositionOverTime(string bodyId)
        {
            var result = new List<(double, Vector3D)>();

            foreach (var frame in _frames)
            {
                var state = frame.BodyStates?.FirstOrDefault(s => s.Id == bodyId);
                if (state != null)
                {
                    result.Add((frame.Time, state.Position));
                }
            }

            return result;
        }

        /// <summary>
        /// Gets velocity over time for a specific body.
        /// </summary>
        public List<(double time, Vector3D velocity)> GetVelocityOverTime(string bodyId)
        {
            var result = new List<(double, Vector3D)>();

            foreach (var frame in _frames)
            {
                var state = frame.BodyStates?.FirstOrDefault(s => s.Id == bodyId);
                if (state != null)
                {
                    result.Add((frame.Time, state.Velocity));
                }
            }

            return result;
        }

        /// <summary>
        /// Calculates the average framerate.
        /// </summary>
        public double CalculateAverageFramerate()
        {
            if (_frames.Count < 2)
                return 0;

            double totalDeltaTime = _frames.Sum(f => f.DeltaTime);
            return _frames.Count / totalDeltaTime;
        }

        #endregion

        #region Export

        /// <summary>
        /// Exports recorded data to CSV format.
        /// </summary>
        public string ExportToCsv()
        {
            var sb = new StringBuilder();

            // Header
            sb.AppendLine("Frame,Time,BodyId,PosX,PosY,PosZ,VelX,VelY,VelZ,KineticEnergy,PotentialEnergy");

            foreach (var frame in _frames)
            {
                if (frame.BodyStates == null)
                    continue;

                foreach (var state in frame.BodyStates)
                {
                    sb.AppendLine(
                        $"{frame.FrameIndex},{frame.Time:F6},{state.Id}," +
                        $"{state.Position.X:F6},{state.Position.Y:F6},{state.Position.Z:F6}," +
                        $"{state.Velocity.X:F6},{state.Velocity.Y:F6},{state.Velocity.Z:F6}," +
                        $"{state.KineticEnergy:F6},{state.PotentialEnergy:F6}"
                    );
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// Exports to CSV file.
        /// </summary>
        public void ExportToCsv(string filePath)
        {
            File.WriteAllText(filePath, ExportToCsv());
        }

        /// <summary>
        /// Exports trajectories to CSV.
        /// </summary>
        public string ExportTrajectoriesToCsv()
        {
            var sb = new StringBuilder();
            sb.AppendLine("BodyId,PointIndex,X,Y,Z");

            foreach (var (bodyId, trajectory) in _trajectories)
            {
                for (int i = 0; i < trajectory.Count; i++)
                {
                    var pos = trajectory[i];
                    sb.AppendLine($"{bodyId},{i},{pos.X:F6},{pos.Y:F6},{pos.Z:F6}");
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// Exports energy data to CSV.
        /// </summary>
        public string ExportEnergyToCsv()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Frame,Time,KineticEnergy,PotentialEnergy,TotalEnergy");

            foreach (var frame in _frames)
            {
                sb.AppendLine(
                    $"{frame.FrameIndex},{frame.Time:F6}," +
                    $"{frame.TotalKineticEnergy:F6},{frame.TotalPotentialEnergy:F6},{frame.TotalEnergy:F6}"
                );
            }

            return sb.ToString();
        }

        /// <summary>
        /// Exports to binary format for efficient storage.
        /// </summary>
        public void ExportToBinary(string filePath)
        {
            using var stream = new FileStream(filePath, FileMode.Create);
            using var writer = new BinaryWriter(stream);

            // Header
            writer.Write("ARTEMIS_REC");
            writer.Write(1); // Version
            writer.Write(_frames.Count);
            writer.Write(_currentTime);

            // Frames
            foreach (var frame in _frames)
            {
                writer.Write(frame.FrameIndex);
                writer.Write(frame.Time);
                writer.Write(frame.DeltaTime);
                writer.Write(frame.BodyStates?.Count ?? 0);

                if (frame.BodyStates != null)
                {
                    foreach (var state in frame.BodyStates)
                    {
                        writer.Write(state.Id);
                        WriteVector3D(writer, state.Position);
                        WriteVector3D(writer, state.Velocity);
                        WriteQuaternion(writer, state.Rotation);
                        WriteVector3D(writer, state.AngularVelocity);
                        writer.Write(state.KineticEnergy);
                        writer.Write(state.PotentialEnergy);
                    }
                }
            }
        }

        /// <summary>
        /// Imports from binary format.
        /// </summary>
        public static SimulationRecorder ImportFromBinary(string filePath)
        {
            var recorder = new SimulationRecorder();

            using var stream = new FileStream(filePath, FileMode.Open);
            using var reader = new BinaryReader(stream);

            // Header
            string magic = reader.ReadString();
            if (magic != "ARTEMIS_REC")
                throw new InvalidDataException("Invalid file format");

            int version = reader.ReadInt32();
            int frameCount = reader.ReadInt32();
            recorder._currentTime = reader.ReadDouble();

            // Frames
            for (int f = 0; f < frameCount; f++)
            {
                var frame = new SimulationFrame
                {
                    FrameIndex = reader.ReadInt32(),
                    Time = reader.ReadDouble(),
                    DeltaTime = reader.ReadDouble(),
                    BodyStates = new List<BodyState>()
                };

                int bodyCount = reader.ReadInt32();
                for (int b = 0; b < bodyCount; b++)
                {
                    frame.BodyStates.Add(new BodyState
                    {
                        Id = reader.ReadString(),
                        Position = ReadVector3D(reader),
                        Velocity = ReadVector3D(reader),
                        Rotation = ReadQuaternion(reader),
                        AngularVelocity = ReadVector3D(reader),
                        KineticEnergy = reader.ReadDouble(),
                        PotentialEnergy = reader.ReadDouble()
                    });
                }

                frame.TotalKineticEnergy = frame.BodyStates.Sum(s => s.KineticEnergy);
                frame.TotalPotentialEnergy = frame.BodyStates.Sum(s => s.PotentialEnergy);
                frame.TotalEnergy = frame.TotalKineticEnergy + frame.TotalPotentialEnergy;

                recorder._frames.Add(frame);
            }

            return recorder;
        }

        private static void WriteVector3D(BinaryWriter writer, Vector3D v)
        {
            writer.Write(v.X);
            writer.Write(v.Y);
            writer.Write(v.Z);
        }

        private static Vector3D ReadVector3D(BinaryReader reader)
        {
            return new Vector3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
        }

        private static void WriteQuaternion(BinaryWriter writer, Quaternion q)
        {
            writer.Write(q.X);
            writer.Write(q.Y);
            writer.Write(q.Z);
            writer.Write(q.W);
        }

        private static Quaternion ReadQuaternion(BinaryReader reader)
        {
            return new Quaternion(
                reader.ReadDouble(),
                reader.ReadDouble(),
                reader.ReadDouble(),
                reader.ReadDouble()
            );
        }

        #endregion
    }

    #region Data Structures

    /// <summary>
    /// Represents a single frame of simulation data.
    /// </summary>
    public class SimulationFrame
    {
        public int FrameIndex { get; set; }
        public double Time { get; set; }
        public double DeltaTime { get; set; }
        public List<BodyState>? BodyStates { get; set; }
        public List<ParticleState>? ParticleStates { get; set; }
        public int ParticleCount { get; set; }
        public double TotalKineticEnergy { get; set; }
        public double TotalPotentialEnergy { get; set; }
        public double TotalEnergy { get; set; }
    }

    /// <summary>
    /// State of a rigid body at a point in time.
    /// </summary>
    public class BodyState
    {
        public string Id { get; set; } = "";
        public Vector3D Position { get; set; }
        public Vector3D Velocity { get; set; }
        public Quaternion Rotation { get; set; }
        public Vector3D AngularVelocity { get; set; }
        public double KineticEnergy { get; set; }
        public double PotentialEnergy { get; set; }
    }

    /// <summary>
    /// State of a particle at a point in time.
    /// </summary>
    public class ParticleState
    {
        public int Index { get; set; }
        public Vector3D Position { get; set; }
        public Vector3D Velocity { get; set; }
        public double Mass { get; set; }
        public double Radius { get; set; }
        public double Lifetime { get; set; }
        public uint Color { get; set; }
    }

    /// <summary>
    /// Statistics for a body over the simulation.
    /// </summary>
    public class BodyStatistics
    {
        public string Id { get; set; } = "";
        public int SampleCount { get; set; }
        public double MaxSpeed { get; set; }
        public double AverageSpeed { get; set; }
        public Vector3D MinPosition { get; set; } = new(double.MaxValue);
        public Vector3D MaxPosition { get; set; } = new(double.MinValue);

        public Vector3D BoundingBoxSize => MaxPosition - MinPosition;
    }

    #endregion
}
