using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;
using Artemis.Destruction;

namespace Artemis.Forces
{
    /// <summary>
    /// G-Force system that tracks acceleration and can destroy objects under extreme G-forces.
    /// </summary>
    public class GForceSystem
    {
        #region Nested Types

        /// <summary>
        /// G-Force limits for different object types.
        /// </summary>
        public class GForceLimits
        {
            /// <summary>Maximum G-force before structural damage begins.</summary>
            public double DamageThreshold { get; set; } = 10.0;

            /// <summary>Maximum G-force before destruction.</summary>
            public double DestructionThreshold { get; set; } = 50.0;

            /// <summary>Direction sensitivity (some objects are weaker in certain directions).</summary>
            public Vector3D DirectionalWeakness { get; set; } = Vector3D.One;

            /// <summary>Time required above threshold before damage occurs.</summary>
            public double DamageDelay { get; set; } = 0.1;
        }

        /// <summary>
        /// Tracking data for a body.
        /// </summary>
        private class BodyGForceData
        {
            public Vector3D PreviousVelocity;
            public double TimeAboveThreshold;
            public double CurrentGForce;
            public double PeakGForce;
            public double AccumulatedDamage;
        }

        #endregion

        #region Fields

        private readonly Dictionary<IPhysicsBody, BodyGForceData> _bodyData;
        private readonly Dictionary<IPhysicsBody, GForceLimits> _bodyLimits;
        private readonly GForceLimits _defaultLimits;
        private readonly List<(IPhysicsBody Body, double GForce, Vector3D Direction)> _destructionQueue;

        // Standard gravity for G calculation (m/s²)
        private const double StandardG = 9.80665;

        #endregion

        #region Properties

        /// <summary>Gets or sets whether the system is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Gets or sets the default G-force limits.</summary>
        public GForceLimits DefaultLimits => _defaultLimits;

        /// <summary>Optional fracture system for destruction effects.</summary>
        public FractureSystem? FractureSystem { get; set; }

        /// <summary>Event raised when a body experiences high G-force.</summary>
        public event Action<IPhysicsBody, double, Vector3D>? OnHighGForce;

        /// <summary>Event raised when a body takes damage from G-force.</summary>
        public event Action<IPhysicsBody, double, double>? OnGForceDamage;

        /// <summary>Event raised when a body is destroyed by G-force.</summary>
        public event Action<IPhysicsBody, double, Vector3D>? OnGForceDestruction;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new G-force system.
        /// </summary>
        public GForceSystem()
        {
            _bodyData = new Dictionary<IPhysicsBody, BodyGForceData>();
            _bodyLimits = new Dictionary<IPhysicsBody, GForceLimits>();
            _defaultLimits = new GForceLimits();
            _destructionQueue = new List<(IPhysicsBody, double, Vector3D)>();
        }

        #endregion

        #region Configuration

        /// <summary>
        /// Registers a body for G-force tracking.
        /// </summary>
        public void RegisterBody(IPhysicsBody body, GForceLimits? limits = null)
        {
            if (!_bodyData.ContainsKey(body))
            {
                _bodyData[body] = new BodyGForceData
                {
                    PreviousVelocity = body.Velocity
                };
            }

            if (limits != null)
            {
                _bodyLimits[body] = limits;
            }
        }

        /// <summary>
        /// Unregisters a body from G-force tracking.
        /// </summary>
        public void UnregisterBody(IPhysicsBody body)
        {
            _bodyData.Remove(body);
            _bodyLimits.Remove(body);
        }

        /// <summary>
        /// Sets custom limits for a body.
        /// </summary>
        public void SetLimits(IPhysicsBody body, GForceLimits limits)
        {
            _bodyLimits[body] = limits;
        }

        #endregion

        #region Update

        /// <summary>
        /// Updates G-force calculations for all tracked bodies.
        /// </summary>
        public void Update(double deltaTime)
        {
            if (!Enabled || deltaTime <= 0) return;

            _destructionQueue.Clear();

            foreach (var kvp in _bodyData)
            {
                var body = kvp.Key;
                var data = kvp.Value;

                if (!body.IsActive) continue;

                // Calculate acceleration (change in velocity / time)
                var acceleration = (body.Velocity - data.PreviousVelocity) / deltaTime;
                double accelerationMagnitude = acceleration.Magnitude;

                // Convert to G-force
                double gForce = accelerationMagnitude / StandardG;
                data.CurrentGForce = gForce;
                data.PeakGForce = Math.Max(data.PeakGForce, gForce);

                // Get limits for this body
                var limits = _bodyLimits.GetValueOrDefault(body, _defaultLimits);

                // Apply directional weakness
                var weakenedAccel = new Vector3D(
                    Math.Abs(acceleration.X) * limits.DirectionalWeakness.X,
                    Math.Abs(acceleration.Y) * limits.DirectionalWeakness.Y,
                    Math.Abs(acceleration.Z) * limits.DirectionalWeakness.Z
                );
                double effectiveGForce = weakenedAccel.Magnitude / StandardG;

                // Check thresholds
                if (effectiveGForce >= limits.DamageThreshold)
                {
                    data.TimeAboveThreshold += deltaTime;
                    OnHighGForce?.Invoke(body, effectiveGForce, acceleration.Normalized);

                    if (data.TimeAboveThreshold >= limits.DamageDelay)
                    {
                        // Calculate damage
                        double damageRate = (effectiveGForce - limits.DamageThreshold) /
                                           (limits.DestructionThreshold - limits.DamageThreshold);
                        damageRate = Math.Clamp(damageRate, 0, 1);
                        double damage = damageRate * deltaTime;
                        data.AccumulatedDamage += damage;

                        OnGForceDamage?.Invoke(body, effectiveGForce, data.AccumulatedDamage);

                        // Check for destruction
                        if (effectiveGForce >= limits.DestructionThreshold || data.AccumulatedDamage >= 1.0)
                        {
                            _destructionQueue.Add((body, effectiveGForce, acceleration.Normalized));
                        }
                    }
                }
                else
                {
                    data.TimeAboveThreshold = 0;
                    // Slowly recover from damage
                    data.AccumulatedDamage = Math.Max(0, data.AccumulatedDamage - deltaTime * 0.1);
                }

                // Store current velocity for next frame
                data.PreviousVelocity = body.Velocity;
            }

            // Process destruction queue
            foreach (var (body, gForce, direction) in _destructionQueue)
            {
                DestroyBody(body, gForce, direction);
            }
        }

        private void DestroyBody(IPhysicsBody body, double gForce, Vector3D direction)
        {
            OnGForceDestruction?.Invoke(body, gForce, direction);

            // If we have a fracture system, use it
            if (FractureSystem != null && body is RigidBody rb)
            {
                var impactPoint = rb.Position + direction * rb.BoundingBox.Size.Magnitude * 0.5;
                var force = gForce * rb.Mass * StandardG * 0.1; // Convert to impulse
                FractureSystem.Fracture(rb, impactPoint, direction, force);
            }

            // Clean up
            UnregisterBody(body);
        }

        #endregion

        #region Queries

        /// <summary>
        /// Gets the current G-force on a body.
        /// </summary>
        public double GetCurrentGForce(IPhysicsBody body)
        {
            return _bodyData.TryGetValue(body, out var data) ? data.CurrentGForce : 0;
        }

        /// <summary>
        /// Gets the peak G-force experienced by a body.
        /// </summary>
        public double GetPeakGForce(IPhysicsBody body)
        {
            return _bodyData.TryGetValue(body, out var data) ? data.PeakGForce : 0;
        }

        /// <summary>
        /// Gets the accumulated damage on a body (0-1).
        /// </summary>
        public double GetAccumulatedDamage(IPhysicsBody body)
        {
            return _bodyData.TryGetValue(body, out var data) ? data.AccumulatedDamage : 0;
        }

        /// <summary>
        /// Resets the peak G-force for a body.
        /// </summary>
        public void ResetPeakGForce(IPhysicsBody body)
        {
            if (_bodyData.TryGetValue(body, out var data))
            {
                data.PeakGForce = 0;
            }
        }

        #endregion

        #region Presets

        /// <summary>Creates limits for human tolerance (6G sustained, 9G brief).</summary>
        public static GForceLimits HumanLimits() => new()
        {
            DamageThreshold = 6.0,
            DestructionThreshold = 9.0,
            DamageDelay = 1.0,
            DirectionalWeakness = new Vector3D(1, 0.7, 1) // Weaker in vertical direction
        };

        /// <summary>Creates limits for trained pilot (9G sustained, 12G brief).</summary>
        public static GForceLimits PilotLimits() => new()
        {
            DamageThreshold = 9.0,
            DestructionThreshold = 12.0,
            DamageDelay = 0.5
        };

        /// <summary>Creates limits for aircraft structure.</summary>
        public static GForceLimits AircraftLimits() => new()
        {
            DamageThreshold = 6.0,
            DestructionThreshold = 9.0,
            DamageDelay = 0.1
        };

        /// <summary>Creates limits for spacecraft.</summary>
        public static GForceLimits SpacecraftLimits() => new()
        {
            DamageThreshold = 4.0,
            DestructionThreshold = 8.0,
            DamageDelay = 0.2
        };

        /// <summary>Creates limits for fragile objects (glass, electronics).</summary>
        public static GForceLimits FragileLimits() => new()
        {
            DamageThreshold = 3.0,
            DestructionThreshold = 10.0,
            DamageDelay = 0.0
        };

        /// <summary>Creates limits for robust objects (tank, armored vehicle).</summary>
        public static GForceLimits RobustLimits() => new()
        {
            DamageThreshold = 50.0,
            DestructionThreshold = 100.0,
            DamageDelay = 0.0
        };

        /// <summary>Creates limits for racing car.</summary>
        public static GForceLimits RaceCarLimits() => new()
        {
            DamageThreshold = 5.0,
            DestructionThreshold = 25.0,
            DamageDelay = 0.05
        };

        #endregion
    }
}
