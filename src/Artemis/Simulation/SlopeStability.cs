using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;
using Artemis.Core;

namespace Artemis.Simulation
{
    /// <summary>
    /// Material properties for slope stability analysis
    /// </summary>
    public class SlopeMaterial
    {
        public string Name { get; set; } = "Generic";
        public float AngleOfRepose { get; set; }         // Degrees - maximum stable slope angle
        public float Cohesion { get; set; }               // Pa - internal cohesion
        public float FrictionAngle { get; set; }          // Degrees - internal friction angle
        public float UnitWeight { get; set; }             // N/m³ - weight per unit volume
        public float WaterContent { get; set; }           // 0-1 - affects stability
        public float WeatheringRate { get; set; }         // Rate of degradation over time
        public float FragmentSizeMin { get; set; } = 0.1f;
        public float FragmentSizeMax { get; set; } = 1.0f;
        public float Bounciness { get; set; } = 0.3f;     // Coefficient of restitution

        // Presets
        public static SlopeMaterial Rock() => new()
        {
            Name = "Rock",
            AngleOfRepose = 45,
            Cohesion = 50000,
            FrictionAngle = 35,
            UnitWeight = 26000,
            FragmentSizeMin = 0.1f,
            FragmentSizeMax = 2.0f,
            Bounciness = 0.4f
        };

        public static SlopeMaterial Gravel() => new()
        {
            Name = "Gravel",
            AngleOfRepose = 35,
            Cohesion = 0,
            FrictionAngle = 38,
            UnitWeight = 18000,
            FragmentSizeMin = 0.02f,
            FragmentSizeMax = 0.1f,
            Bounciness = 0.2f
        };

        public static SlopeMaterial Sand() => new()
        {
            Name = "Sand",
            AngleOfRepose = 34,
            Cohesion = 0,
            FrictionAngle = 32,
            UnitWeight = 16000,
            FragmentSizeMin = 0.001f,
            FragmentSizeMax = 0.002f,
            Bounciness = 0.1f
        };

        public static SlopeMaterial Clay() => new()
        {
            Name = "Clay",
            AngleOfRepose = 25,
            Cohesion = 20000,
            FrictionAngle = 20,
            UnitWeight = 19000,
            WaterContent = 0.3f,
            FragmentSizeMin = 0.05f,
            FragmentSizeMax = 0.3f,
            Bounciness = 0.05f
        };

        public static SlopeMaterial Soil() => new()
        {
            Name = "Soil",
            AngleOfRepose = 30,
            Cohesion = 5000,
            FrictionAngle = 28,
            UnitWeight = 17000,
            WaterContent = 0.15f,
            FragmentSizeMin = 0.01f,
            FragmentSizeMax = 0.1f,
            Bounciness = 0.1f
        };

        public static SlopeMaterial Ice() => new()
        {
            Name = "Ice",
            AngleOfRepose = 35,
            Cohesion = 100000,
            FrictionAngle = 15,
            UnitWeight = 9000,
            FragmentSizeMin = 0.05f,
            FragmentSizeMax = 0.5f,
            Bounciness = 0.5f
        };

        public static SlopeMaterial Snow() => new()
        {
            Name = "Snow",
            AngleOfRepose = 38,
            Cohesion = 1000,
            FrictionAngle = 25,
            UnitWeight = 3000,
            WaterContent = 0.1f,
            FragmentSizeMin = 0.005f,
            FragmentSizeMax = 0.05f,
            Bounciness = 0.05f
        };
    }

    /// <summary>
    /// Represents a point on a slope that can become unstable
    /// </summary>
    public class SlopeCell
    {
        public Vector3 Position { get; set; }
        public Vector3 Normal { get; set; } = Vector3.UnitY;
        public float Height { get; set; }
        public float Mass { get; set; }
        public SlopeMaterial Material { get; set; } = SlopeMaterial.Rock();
        public float StabilityFactor { get; set; } = 1.0f;  // >1 = stable, <1 = unstable
        public float Saturation { get; set; }  // Water saturation 0-1
        public bool IsUnstable => StabilityFactor < 1.0f;
        public bool HasCollapsed { get; set; }
        public float WeatheringProgress { get; set; }  // 0-1, increases over time
    }

    /// <summary>
    /// Rockfall fragment during simulation
    /// </summary>
    public class RockFragment
    {
        public Vector3 Position { get; set; }
        public Vector3 Velocity { get; set; }
        public Vector3 AngularVelocity { get; set; }
        public float Mass { get; set; }
        public float Radius { get; set; }
        public float Bounciness { get; set; }
        public bool IsActive { get; set; } = true;
        public float Energy => 0.5f * Mass * Velocity.LengthSquared();
        public int BounceCount { get; set; }
        public float TravelDistance { get; set; }
        public List<Vector3> Trajectory { get; } = new();
    }

    /// <summary>
    /// Result of a rockfall simulation
    /// </summary>
    public class RockfallResult
    {
        public List<RockFragment> Fragments { get; } = new();
        public List<Vector3> ImpactPoints { get; } = new();
        public float MaxRunoutDistance { get; set; }
        public float MaxEnergy { get; set; }
        public float TotalMassReleased { get; set; }
        public float SimulationTime { get; set; }
        public bool IsComplete { get; set; }
    }

    /// <summary>
    /// Precomputed rockfall scenario
    /// </summary>
    public class PrecomputedRockfall
    {
        public Vector3 SourcePosition { get; set; }
        public float SourceMass { get; set; }
        public SlopeMaterial Material { get; set; } = SlopeMaterial.Rock();
        public List<RockfallFrame> Frames { get; } = new();
        public float TotalDuration { get; set; }
        public RockfallResult FinalResult { get; set; } = new();

        /// <summary>
        /// Play precomputed simulation forward
        /// </summary>
        public RockfallFrame? GetFrame(float time)
        {
            if (Frames.Count == 0) return null;
            int index = (int)(time / TotalDuration * Frames.Count);
            index = Math.Clamp(index, 0, Frames.Count - 1);
            return Frames[index];
        }

        /// <summary>
        /// Play precomputed simulation in reverse (time reversal)
        /// </summary>
        public RockfallFrame? GetFrameReverse(float time)
        {
            if (Frames.Count == 0) return null;
            int index = Frames.Count - 1 - (int)(time / TotalDuration * Frames.Count);
            index = Math.Clamp(index, 0, Frames.Count - 1);

            // Reverse velocities for backwards playback
            var frame = Frames[index];
            var reversed = new RockfallFrame { Time = time };
            foreach (var frag in frame.FragmentStates)
            {
                reversed.FragmentStates.Add(new FragmentState
                {
                    Position = frag.Position,
                    Velocity = -frag.Velocity,  // Reverse velocity
                    Mass = frag.Mass,
                    Radius = frag.Radius,
                    IsActive = frag.IsActive
                });
            }
            return reversed;
        }
    }

    public class RockfallFrame
    {
        public float Time { get; set; }
        public List<FragmentState> FragmentStates { get; } = new();
    }

    public class FragmentState
    {
        public Vector3 Position { get; set; }
        public Vector3 Velocity { get; set; }
        public float Mass { get; set; }
        public float Radius { get; set; }
        public bool IsActive { get; set; }
    }

    /// <summary>
    /// Slope stability analysis and rockfall simulation
    /// </summary>
    public class SlopeStabilitySystem
    {
        private readonly List<SlopeCell> _cells = new();
        private readonly List<RockFragment> _activeFragments = new();
        private readonly List<PrecomputedRockfall> _precomputed = new();
        private Func<Vector3, float>? _terrainHeightFunc;
        private Func<Vector3, Vector3>? _terrainNormalFunc;

        public Vector3 Gravity { get; set; } = new(0, -9.81f, 0);
        public float TimeScale { get; set; } = 1.0f;
        public int MaxFragments { get; set; } = 1000;
        public float MinFragmentEnergy { get; set; } = 0.1f;  // Stop when energy drops below this
        public int MaxBounces { get; set; } = 50;
        public bool RecordTrajectories { get; set; } = true;

        public IReadOnlyList<SlopeCell> Cells => _cells;
        public IReadOnlyList<RockFragment> ActiveFragments => _activeFragments;
        public IReadOnlyList<PrecomputedRockfall> PrecomputedRockfalls => _precomputed;

        /// <summary>
        /// Set terrain height function for collision detection
        /// </summary>
        public void SetTerrainFunction(Func<Vector3, float> heightFunc, Func<Vector3, Vector3>? normalFunc = null)
        {
            _terrainHeightFunc = heightFunc;
            _terrainNormalFunc = normalFunc ?? (pos =>
            {
                // Estimate normal from height gradient
                float eps = 0.1f;
                float hx = heightFunc(pos + Vector3.UnitX * eps) - heightFunc(pos - Vector3.UnitX * eps);
                float hz = heightFunc(pos + Vector3.UnitZ * eps) - heightFunc(pos - Vector3.UnitZ * eps);
                return Vector3.Normalize(new Vector3(-hx, 2 * eps, -hz));
            });
        }

        /// <summary>
        /// Add a slope cell for stability analysis
        /// </summary>
        public void AddCell(SlopeCell cell)
        {
            _cells.Add(cell);
        }

        /// <summary>
        /// Analyze stability of a single slope cell using Factor of Safety
        /// </summary>
        public float AnalyzeStability(SlopeCell cell)
        {
            // Calculate slope angle from normal
            float slopeAngle = MathF.Acos(Vector3.Dot(cell.Normal, Vector3.UnitY)) * 180f / MathF.PI;

            // Mohr-Coulomb failure criterion
            // Factor of Safety = (c + σn * tan(φ)) / τ
            // Where:
            // c = cohesion
            // σn = normal stress
            // φ = friction angle
            // τ = shear stress

            float phi = cell.Material.FrictionAngle * MathF.PI / 180f;
            float slopeRad = slopeAngle * MathF.PI / 180f;

            // Normal and shear stresses on failure plane
            float normalStress = cell.Material.UnitWeight * cell.Height * MathF.Cos(slopeRad) * MathF.Cos(slopeRad);
            float shearStress = cell.Material.UnitWeight * cell.Height * MathF.Sin(slopeRad) * MathF.Cos(slopeRad);

            // Water pressure reduces effective stress
            float waterPressure = cell.Saturation * 9810 * cell.Height; // γw * h
            float effectiveNormalStress = MathF.Max(0, normalStress - waterPressure);

            // Shear strength
            float shearStrength = cell.Material.Cohesion + effectiveNormalStress * MathF.Tan(phi);

            // Factor of Safety
            float factorOfSafety = shearStress > 0 ? shearStrength / shearStress : float.MaxValue;

            // Apply weathering reduction
            factorOfSafety *= (1 - cell.WeatheringProgress * 0.5f);

            cell.StabilityFactor = factorOfSafety;
            return factorOfSafety;
        }

        /// <summary>
        /// Analyze all cells
        /// </summary>
        public void AnalyzeAllCells()
        {
            Parallel.ForEach(_cells, cell => AnalyzeStability(cell));
        }

        /// <summary>
        /// Find unstable cells
        /// </summary>
        public List<SlopeCell> FindUnstableCells(float threshold = 1.0f)
        {
            return _cells.Where(c => c.StabilityFactor < threshold && !c.HasCollapsed).ToList();
        }

        /// <summary>
        /// Trigger collapse at a cell, generating rock fragments
        /// </summary>
        public List<RockFragment> TriggerCollapse(SlopeCell cell, int fragmentCount = 10)
        {
            if (cell.HasCollapsed) return new List<RockFragment>();

            cell.HasCollapsed = true;
            var fragments = new List<RockFragment>();
            var random = new Random();

            float totalMass = cell.Mass;
            float remainingMass = totalMass;

            for (int i = 0; i < fragmentCount && remainingMass > 0; i++)
            {
                float fragMass = totalMass / fragmentCount * (0.5f + (float)random.NextDouble());
                fragMass = MathF.Min(fragMass, remainingMass);
                remainingMass -= fragMass;

                float radius = cell.Material.FragmentSizeMin +
                    (float)random.NextDouble() * (cell.Material.FragmentSizeMax - cell.Material.FragmentSizeMin);

                var fragment = new RockFragment
                {
                    Position = cell.Position + new Vector3(
                        ((float)random.NextDouble() - 0.5f) * radius * 2,
                        radius,
                        ((float)random.NextDouble() - 0.5f) * radius * 2
                    ),
                    Velocity = GetInitialVelocity(cell, random),
                    Mass = fragMass,
                    Radius = radius,
                    Bounciness = cell.Material.Bounciness
                };

                fragments.Add(fragment);
                if (_activeFragments.Count < MaxFragments)
                {
                    _activeFragments.Add(fragment);
                }
            }

            return fragments;
        }

        private Vector3 GetInitialVelocity(SlopeCell cell, Random random)
        {
            // Initial velocity is down-slope plus some randomness
            Vector3 downSlope = Vector3.Cross(Vector3.Cross(cell.Normal, Vector3.UnitY), cell.Normal);
            if (downSlope.LengthSquared() < 0.001f)
            {
                downSlope = new Vector3((float)random.NextDouble() - 0.5f, 0, (float)random.NextDouble() - 0.5f);
            }
            downSlope = Vector3.Normalize(downSlope);

            float speed = 1f + (float)random.NextDouble() * 2f;
            Vector3 randomComponent = new Vector3(
                ((float)random.NextDouble() - 0.5f) * 0.5f,
                (float)random.NextDouble() * 0.5f,
                ((float)random.NextDouble() - 0.5f) * 0.5f
            );

            return downSlope * speed + randomComponent;
        }

        /// <summary>
        /// Update real-time rockfall simulation
        /// </summary>
        public void Update(float deltaTime)
        {
            deltaTime *= TimeScale;

            for (int i = _activeFragments.Count - 1; i >= 0; i--)
            {
                var frag = _activeFragments[i];
                if (!frag.IsActive) continue;

                // Apply gravity
                frag.Velocity += Gravity * deltaTime;

                // Update position
                Vector3 newPos = frag.Position + frag.Velocity * deltaTime;

                // Check terrain collision
                if (_terrainHeightFunc != null)
                {
                    float terrainHeight = _terrainHeightFunc(newPos);

                    if (newPos.Y - frag.Radius < terrainHeight)
                    {
                        // Collision with terrain
                        Vector3 normal = _terrainNormalFunc?.Invoke(newPos) ?? Vector3.UnitY;

                        // Reflect velocity with energy loss
                        float vn = Vector3.Dot(frag.Velocity, normal);
                        if (vn < 0)  // Moving into surface
                        {
                            frag.Velocity -= (1 + frag.Bounciness) * vn * normal;

                            // Friction loss
                            Vector3 tangent = frag.Velocity - Vector3.Dot(frag.Velocity, normal) * normal;
                            frag.Velocity = tangent * 0.8f + Vector3.Dot(frag.Velocity, normal) * normal;

                            frag.BounceCount++;

                            // Add angular velocity from impact
                            frag.AngularVelocity += Vector3.Cross(normal, frag.Velocity) * 0.1f;
                        }

                        // Place on surface
                        newPos.Y = terrainHeight + frag.Radius;
                    }
                }

                frag.Position = newPos;
                frag.TravelDistance += (frag.Velocity * deltaTime).Length();

                if (RecordTrajectories && frag.Trajectory.Count < 1000)
                {
                    frag.Trajectory.Add(frag.Position);
                }

                // Check if fragment should stop
                if (frag.Energy < MinFragmentEnergy || frag.BounceCount > MaxBounces)
                {
                    frag.IsActive = false;
                }
            }

            // Remove inactive fragments if too many
            if (_activeFragments.Count(f => !f.IsActive) > MaxFragments / 2)
            {
                _activeFragments.RemoveAll(f => !f.IsActive);
            }
        }

        /// <summary>
        /// Precompute a rockfall scenario
        /// </summary>
        public PrecomputedRockfall PrecomputeRockfall(
            Vector3 sourcePosition,
            float sourceMass,
            SlopeMaterial material,
            float duration,
            float timeStep = 0.02f,
            int fragmentCount = 20)
        {
            var precomputed = new PrecomputedRockfall
            {
                SourcePosition = sourcePosition,
                SourceMass = sourceMass,
                Material = material,
                TotalDuration = duration
            };

            // Create source cell
            var sourceCell = new SlopeCell
            {
                Position = sourcePosition,
                Mass = sourceMass,
                Material = material,
                Normal = _terrainNormalFunc?.Invoke(sourcePosition) ?? Vector3.UnitY
            };

            // Generate initial fragments
            var fragments = new List<RockFragment>();
            var random = new Random(42);  // Fixed seed for reproducibility

            float remainingMass = sourceMass;
            for (int i = 0; i < fragmentCount && remainingMass > 0; i++)
            {
                float fragMass = sourceMass / fragmentCount * (0.5f + (float)random.NextDouble());
                fragMass = MathF.Min(fragMass, remainingMass);
                remainingMass -= fragMass;

                float radius = material.FragmentSizeMin +
                    (float)random.NextDouble() * (material.FragmentSizeMax - material.FragmentSizeMin);

                fragments.Add(new RockFragment
                {
                    Position = sourcePosition + new Vector3(
                        ((float)random.NextDouble() - 0.5f) * radius * 2,
                        radius,
                        ((float)random.NextDouble() - 0.5f) * radius * 2
                    ),
                    Velocity = GetInitialVelocity(sourceCell, random),
                    Mass = fragMass,
                    Radius = radius,
                    Bounciness = material.Bounciness,
                    IsActive = true
                });

                precomputed.FinalResult.TotalMassReleased += fragMass;
            }

            // Simulate
            float time = 0;
            while (time < duration)
            {
                // Record frame
                var frame = new RockfallFrame { Time = time };
                foreach (var frag in fragments)
                {
                    frame.FragmentStates.Add(new FragmentState
                    {
                        Position = frag.Position,
                        Velocity = frag.Velocity,
                        Mass = frag.Mass,
                        Radius = frag.Radius,
                        IsActive = frag.IsActive
                    });
                }
                precomputed.Frames.Add(frame);

                // Update fragments
                foreach (var frag in fragments)
                {
                    if (!frag.IsActive) continue;

                    // Apply gravity
                    frag.Velocity += Gravity * timeStep;
                    Vector3 newPos = frag.Position + frag.Velocity * timeStep;

                    // Terrain collision
                    if (_terrainHeightFunc != null)
                    {
                        float terrainHeight = _terrainHeightFunc(newPos);
                        if (newPos.Y - frag.Radius < terrainHeight)
                        {
                            Vector3 normal = _terrainNormalFunc?.Invoke(newPos) ?? Vector3.UnitY;
                            float vn = Vector3.Dot(frag.Velocity, normal);
                            if (vn < 0)
                            {
                                frag.Velocity -= (1 + frag.Bounciness) * vn * normal;
                                Vector3 tangent = frag.Velocity - Vector3.Dot(frag.Velocity, normal) * normal;
                                frag.Velocity = tangent * 0.8f + Vector3.Dot(frag.Velocity, normal) * normal;
                                frag.BounceCount++;

                                precomputed.FinalResult.ImpactPoints.Add(newPos);
                            }
                            newPos.Y = terrainHeight + frag.Radius;
                        }
                    }

                    frag.Position = newPos;
                    frag.TravelDistance += (frag.Velocity * timeStep).Length();

                    if (frag.Energy < MinFragmentEnergy || frag.BounceCount > MaxBounces)
                    {
                        frag.IsActive = false;
                    }

                    // Track stats
                    if (frag.TravelDistance > precomputed.FinalResult.MaxRunoutDistance)
                    {
                        precomputed.FinalResult.MaxRunoutDistance = frag.TravelDistance;
                    }
                    if (frag.Energy > precomputed.FinalResult.MaxEnergy)
                    {
                        precomputed.FinalResult.MaxEnergy = frag.Energy;
                    }
                }

                time += timeStep;
            }

            // Final state
            foreach (var frag in fragments)
            {
                precomputed.FinalResult.Fragments.Add(frag);
            }
            precomputed.FinalResult.SimulationTime = duration;
            precomputed.FinalResult.IsComplete = true;

            _precomputed.Add(precomputed);
            return precomputed;
        }

        /// <summary>
        /// Apply weathering to all cells over time
        /// </summary>
        public void ApplyWeathering(float deltaTime)
        {
            foreach (var cell in _cells)
            {
                cell.WeatheringProgress += cell.Material.WeatheringRate * deltaTime;
                cell.WeatheringProgress = MathF.Min(1, cell.WeatheringProgress);
            }
        }

        /// <summary>
        /// Apply rainfall/saturation to slope
        /// </summary>
        public void ApplyRainfall(float intensity, float deltaTime)
        {
            foreach (var cell in _cells)
            {
                cell.Saturation = MathF.Min(1, cell.Saturation + intensity * deltaTime);
            }
        }

        /// <summary>
        /// Apply drainage to slope
        /// </summary>
        public void ApplyDrainage(float rate, float deltaTime)
        {
            foreach (var cell in _cells)
            {
                cell.Saturation = MathF.Max(0, cell.Saturation - rate * deltaTime);
            }
        }

        /// <summary>
        /// Check for chain reaction collapses
        /// </summary>
        public List<SlopeCell> PropagateCollapse(SlopeCell initialCell, float propagationRadius = 5f)
        {
            var collapsed = new List<SlopeCell> { initialCell };
            var queue = new Queue<SlopeCell>();
            queue.Enqueue(initialCell);

            while (queue.Count > 0)
            {
                var current = queue.Dequeue();

                foreach (var cell in _cells)
                {
                    if (cell.HasCollapsed) continue;
                    if ((cell.Position - current.Position).Length() > propagationRadius) continue;

                    // Nearby collapse reduces stability
                    cell.StabilityFactor *= 0.7f;

                    if (cell.StabilityFactor < 1.0f)
                    {
                        cell.HasCollapsed = true;
                        collapsed.Add(cell);
                        queue.Enqueue(cell);
                    }
                }
            }

            return collapsed;
        }

        /// <summary>
        /// Get runout zone boundary for hazard mapping
        /// </summary>
        public List<Vector3> GetRunoutZoneBoundary(PrecomputedRockfall rockfall, float percentile = 0.95f)
        {
            var positions = rockfall.FinalResult.Fragments
                .Select(f => f.Position)
                .OrderBy(p => (p - rockfall.SourcePosition).Length())
                .ToList();

            int boundaryIndex = (int)(positions.Count * percentile);
            return positions.Take(boundaryIndex + 1).ToList();
        }

        /// <summary>
        /// Clear all simulation data
        /// </summary>
        public void Clear()
        {
            _cells.Clear();
            _activeFragments.Clear();
        }
    }
}
