using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Text.Json;
using System.Text.Json.Serialization;
using Artemis.Core;
using Artemis.Particles;

namespace Artemis.Simulation
{
    /// <summary>
    /// Serializable physics state snapshot
    /// </summary>
    public class PhysicsState
    {
        public float Time { get; set; }
        public string Name { get; set; } = "State";
        public List<RigidBodyState> Bodies { get; set; } = new();
        public List<ParticleSnapshot> Particles { get; set; } = new();
        public Dictionary<string, object> CustomData { get; set; } = new();
        public StateMetadata Metadata { get; set; } = new();
    }

    public class StateMetadata
    {
        public DateTime CreatedAt { get; set; } = DateTime.UtcNow;
        public string Version { get; set; } = "1.0";
        public string Description { get; set; } = "";
        public bool IsInitialState { get; set; }
        public bool IsFinalState { get; set; }
        public float SimulationDuration { get; set; }
        public int FrameCount { get; set; }
    }

    public class RigidBodyState
    {
        public int Id { get; set; }
        public string Name { get; set; } = "";
        public Vector3Serializable Position { get; set; }
        public Vector3Serializable Velocity { get; set; }
        public Vector3Serializable AngularVelocity { get; set; }
        public QuaternionSerializable Rotation { get; set; }
        public double Mass { get; set; }
        public bool IsStatic { get; set; }
        public bool IsActive { get; set; } = true;
        public Dictionary<string, float> Properties { get; set; } = new();
    }

    public class ParticleSnapshot
    {
        public Vector3Serializable Position { get; set; }
        public Vector3Serializable Velocity { get; set; }
        public double Mass { get; set; }
        public double Lifetime { get; set; }
        public bool IsAlive { get; set; }
        public int Type { get; set; }
    }

    // JSON-serializable Vector3D
    public struct Vector3Serializable
    {
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }

        public Vector3Serializable(double x, double y, double z) { X = x; Y = y; Z = z; }
        public Vector3Serializable(Vector3D v) { X = v.X; Y = v.Y; Z = v.Z; }
        public Vector3D ToVector3D() => new(X, Y, Z);
        public static implicit operator Vector3Serializable(Vector3D v) => new(v);
        public static implicit operator Vector3D(Vector3Serializable v) => v.ToVector3D();
    }

    // JSON-serializable Quaternion
    public struct QuaternionSerializable
    {
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double W { get; set; }

        public QuaternionSerializable(double x, double y, double z, double w) { X = x; Y = y; Z = z; W = w; }
        public QuaternionSerializable(Quaternion q) { X = q.X; Y = q.Y; Z = q.Z; W = q.W; }
        public Quaternion ToQuaternion() => new(X, Y, Z, W);
        public static implicit operator QuaternionSerializable(Quaternion q) => new(q);
        public static implicit operator Quaternion(QuaternionSerializable q) => q.ToQuaternion();
    }

    /// <summary>
    /// Manages state capture, transfer, and time reversal
    /// </summary>
    public class StateTransferSystem
    {
        private readonly List<PhysicsState> _stateHistory = new();
        private readonly Dictionary<string, PhysicsState> _namedStates = new();
        private int _maxHistorySize = 1000;
        private int _currentStateIndex = -1;

        public int MaxHistorySize
        {
            get => _maxHistorySize;
            set => _maxHistorySize = Math.Max(1, value);
        }

        public int StateCount => _stateHistory.Count;
        public bool CanUndo => _currentStateIndex > 0;
        public bool CanRedo => _currentStateIndex < _stateHistory.Count - 1;
        public PhysicsState? CurrentState => _currentStateIndex >= 0 && _currentStateIndex < _stateHistory.Count
            ? _stateHistory[_currentStateIndex]
            : null;

        /// <summary>
        /// Capture current state of rigid bodies
        /// </summary>
        public PhysicsState CaptureState(IEnumerable<RigidBody> bodies, float time, string name = "")
        {
            var state = new PhysicsState
            {
                Time = time,
                Name = string.IsNullOrEmpty(name) ? $"State_{time:F2}" : name
            };

            int id = 0;
            foreach (var body in bodies)
            {
                state.Bodies.Add(new RigidBodyState
                {
                    Id = id++,
                    Position = body.Position,
                    Velocity = body.Velocity,
                    AngularVelocity = body.AngularVelocity,
                    Rotation = body.Rotation,
                    Mass = body.Mass,
                    IsStatic = body.BodyType == BodyType.Static,
                    IsActive = body.IsActive
                });
            }

            return state;
        }

        /// <summary>
        /// Capture current state of particles
        /// </summary>
        public PhysicsState CaptureState(Span<Particle> particles, float time, string name = "")
        {
            var state = new PhysicsState
            {
                Time = time,
                Name = string.IsNullOrEmpty(name) ? $"State_{time:F2}" : name
            };

            foreach (ref readonly var particle in particles)
            {
                if (!particle.IsAlive) continue;

                state.Particles.Add(new ParticleSnapshot
                {
                    Position = particle.Position,
                    Velocity = particle.Velocity,
                    Mass = particle.Mass,
                    Lifetime = particle.Lifetime,
                    IsAlive = true
                });
            }

            return state;
        }

        /// <summary>
        /// Capture combined state
        /// </summary>
        public PhysicsState CaptureState(IEnumerable<RigidBody> bodies, Span<Particle> particles, float time, string name = "")
        {
            var state = CaptureState(bodies, time, name);
            foreach (ref readonly var particle in particles)
            {
                if (!particle.IsAlive) continue;

                state.Particles.Add(new ParticleSnapshot
                {
                    Position = particle.Position,
                    Velocity = particle.Velocity,
                    Mass = particle.Mass,
                    Lifetime = particle.Lifetime,
                    IsAlive = true
                });
            }
            return state;
        }

        /// <summary>
        /// Apply state to rigid bodies
        /// </summary>
        public void ApplyState(PhysicsState state, List<RigidBody> bodies)
        {
            // Match by index (simple case) or create new bodies
            while (bodies.Count < state.Bodies.Count)
            {
                bodies.Add(new RigidBody());
            }

            for (int i = 0; i < state.Bodies.Count; i++)
            {
                var bs = state.Bodies[i];
                var body = bodies[i];

                body.Position = bs.Position;
                body.Velocity = bs.Velocity;
                body.AngularVelocity = bs.AngularVelocity;
                body.Rotation = bs.Rotation;
                body.Mass = bs.Mass;
                body.BodyType = bs.IsStatic ? BodyType.Static : BodyType.Dynamic;
                body.IsActive = bs.IsActive;
            }
        }

        /// <summary>
        /// Apply state to particles
        /// </summary>
        public void ApplyState(PhysicsState state, Span<Particle> particles)
        {
            int count = Math.Min(state.Particles.Count, particles.Length);
            for (int i = 0; i < count; i++)
            {
                var ps = state.Particles[i];
                particles[i] = new Particle
                {
                    Position = ps.Position,
                    Velocity = ps.Velocity,
                    Mass = ps.Mass,
                    Lifetime = ps.Lifetime,
                    InitialLifetime = ps.Lifetime,
                    IsAlive = ps.IsAlive
                };
            }
        }

        /// <summary>
        /// Create time-reversed state (velocities negated)
        /// </summary>
        public PhysicsState CreateReversedState(PhysicsState state)
        {
            var reversed = new PhysicsState
            {
                Time = state.Time,
                Name = $"{state.Name}_reversed",
                Metadata = new StateMetadata
                {
                    Description = $"Time-reversed: {state.Metadata.Description}",
                    IsInitialState = state.Metadata.IsFinalState,
                    IsFinalState = state.Metadata.IsInitialState
                }
            };

            foreach (var bs in state.Bodies)
            {
                reversed.Bodies.Add(new RigidBodyState
                {
                    Id = bs.Id,
                    Name = bs.Name,
                    Position = bs.Position,
                    Velocity = new Vector3Serializable(-bs.Velocity.X, -bs.Velocity.Y, -bs.Velocity.Z),
                    AngularVelocity = new Vector3Serializable(-bs.AngularVelocity.X, -bs.AngularVelocity.Y, -bs.AngularVelocity.Z),
                    Rotation = bs.Rotation,
                    Mass = bs.Mass,
                    IsStatic = bs.IsStatic,
                    IsActive = bs.IsActive,
                    Properties = bs.Properties
                });
            }

            foreach (var ps in state.Particles)
            {
                reversed.Particles.Add(new ParticleSnapshot
                {
                    Position = ps.Position,
                    Velocity = new Vector3Serializable(-ps.Velocity.X, -ps.Velocity.Y, -ps.Velocity.Z),
                    Mass = ps.Mass,
                    Lifetime = ps.Lifetime,
                    IsAlive = ps.IsAlive,
                    Type = ps.Type
                });
            }

            return reversed;
        }

        /// <summary>
        /// Use final state as initial state (with optional velocity reversal)
        /// </summary>
        public PhysicsState FinalToInitial(PhysicsState finalState, bool reverseVelocities = false)
        {
            var initial = reverseVelocities ? CreateReversedState(finalState) : CloneState(finalState);
            initial.Metadata.IsInitialState = true;
            initial.Metadata.IsFinalState = false;
            initial.Time = 0;
            initial.Name = $"{finalState.Name}_as_initial";
            return initial;
        }

        /// <summary>
        /// Clone a state
        /// </summary>
        public PhysicsState CloneState(PhysicsState state)
        {
            var json = JsonSerializer.Serialize(state);
            return JsonSerializer.Deserialize<PhysicsState>(json) ?? new PhysicsState();
        }

        /// <summary>
        /// Add state to history (for undo/redo)
        /// </summary>
        public void PushState(PhysicsState state)
        {
            // Remove any states after current position (discard redo history)
            if (_currentStateIndex < _stateHistory.Count - 1)
            {
                _stateHistory.RemoveRange(_currentStateIndex + 1, _stateHistory.Count - _currentStateIndex - 1);
            }

            _stateHistory.Add(state);
            _currentStateIndex = _stateHistory.Count - 1;

            // Trim if exceeds max size
            if (_stateHistory.Count > _maxHistorySize)
            {
                _stateHistory.RemoveAt(0);
                _currentStateIndex--;
            }
        }

        /// <summary>
        /// Get previous state (undo)
        /// </summary>
        public PhysicsState? Undo()
        {
            if (!CanUndo) return null;
            _currentStateIndex--;
            return _stateHistory[_currentStateIndex];
        }

        /// <summary>
        /// Get next state (redo)
        /// </summary>
        public PhysicsState? Redo()
        {
            if (!CanRedo) return null;
            _currentStateIndex++;
            return _stateHistory[_currentStateIndex];
        }

        /// <summary>
        /// Store state by name
        /// </summary>
        public void SaveNamedState(string name, PhysicsState state)
        {
            _namedStates[name] = state;
        }

        /// <summary>
        /// Retrieve state by name
        /// </summary>
        public PhysicsState? LoadNamedState(string name)
        {
            return _namedStates.TryGetValue(name, out var state) ? state : null;
        }

        /// <summary>
        /// Interpolate between two states
        /// </summary>
        public PhysicsState Interpolate(PhysicsState stateA, PhysicsState stateB, float t)
        {
            t = Math.Clamp(t, 0, 1);

            var result = new PhysicsState
            {
                Time = stateA.Time + (stateB.Time - stateA.Time) * t,
                Name = $"Interpolated_{t:F2}"
            };

            int bodyCount = Math.Min(stateA.Bodies.Count, stateB.Bodies.Count);
            for (int i = 0; i < bodyCount; i++)
            {
                var a = stateA.Bodies[i];
                var b = stateB.Bodies[i];

                result.Bodies.Add(new RigidBodyState
                {
                    Id = a.Id,
                    Name = a.Name,
                    Position = Vector3D.Lerp(a.Position, b.Position, t),
                    Velocity = Vector3D.Lerp(a.Velocity, b.Velocity, t),
                    AngularVelocity = Vector3D.Lerp(a.AngularVelocity, b.AngularVelocity, t),
                    Rotation = Quaternion.Slerp(a.Rotation, b.Rotation, t),
                    Mass = a.Mass + (b.Mass - a.Mass) * t,
                    IsStatic = a.IsStatic,
                    IsActive = a.IsActive || b.IsActive
                });
            }

            int particleCount = Math.Min(stateA.Particles.Count, stateB.Particles.Count);
            for (int i = 0; i < particleCount; i++)
            {
                var a = stateA.Particles[i];
                var b = stateB.Particles[i];

                result.Particles.Add(new ParticleSnapshot
                {
                    Position = Vector3D.Lerp(a.Position, b.Position, t),
                    Velocity = Vector3D.Lerp(a.Velocity, b.Velocity, t),
                    Mass = a.Mass + (b.Mass - a.Mass) * t,
                    Lifetime = a.Lifetime + (b.Lifetime - a.Lifetime) * t,
                    IsAlive = a.IsAlive || b.IsAlive
                });
            }

            return result;
        }

        /// <summary>
        /// Save state to file
        /// </summary>
        public void SaveToFile(PhysicsState state, string filePath, bool compress = true)
        {
            var json = JsonSerializer.Serialize(state, new JsonSerializerOptions { WriteIndented = !compress });

            if (compress)
            {
                using var fileStream = File.Create(filePath);
                using var gzipStream = new GZipStream(fileStream, CompressionLevel.Optimal);
                using var writer = new StreamWriter(gzipStream);
                writer.Write(json);
            }
            else
            {
                File.WriteAllText(filePath, json);
            }
        }

        /// <summary>
        /// Load state from file
        /// </summary>
        public PhysicsState? LoadFromFile(string filePath)
        {
            string json;

            try
            {
                // Try compressed first
                using var fileStream = File.OpenRead(filePath);
                using var gzipStream = new GZipStream(fileStream, CompressionMode.Decompress);
                using var reader = new StreamReader(gzipStream);
                json = reader.ReadToEnd();
            }
            catch
            {
                // Fall back to uncompressed
                json = File.ReadAllText(filePath);
            }

            return JsonSerializer.Deserialize<PhysicsState>(json);
        }

        /// <summary>
        /// Create a trajectory from recorded states
        /// </summary>
        public List<PhysicsState> CreateTrajectory(PhysicsState initial, PhysicsState final, int frameCount)
        {
            var trajectory = new List<PhysicsState>();

            for (int i = 0; i <= frameCount; i++)
            {
                float t = (float)i / frameCount;
                trajectory.Add(Interpolate(initial, final, t));
            }

            return trajectory;
        }

        /// <summary>
        /// Play trajectory forward
        /// </summary>
        public IEnumerable<PhysicsState> PlayForward(List<PhysicsState> trajectory)
        {
            foreach (var state in trajectory)
            {
                yield return state;
            }
        }

        /// <summary>
        /// Play trajectory in reverse
        /// </summary>
        public IEnumerable<PhysicsState> PlayReverse(List<PhysicsState> trajectory)
        {
            for (int i = trajectory.Count - 1; i >= 0; i--)
            {
                yield return CreateReversedState(trajectory[i]);
            }
        }

        /// <summary>
        /// Clear all history
        /// </summary>
        public void ClearHistory()
        {
            _stateHistory.Clear();
            _currentStateIndex = -1;
        }

        /// <summary>
        /// Export all named states
        /// </summary>
        public Dictionary<string, PhysicsState> ExportNamedStates()
        {
            return new Dictionary<string, PhysicsState>(_namedStates);
        }

        /// <summary>
        /// Import named states
        /// </summary>
        public void ImportNamedStates(Dictionary<string, PhysicsState> states)
        {
            foreach (var kvp in states)
            {
                _namedStates[kvp.Key] = kvp.Value;
            }
        }
    }

}
