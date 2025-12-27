using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Forces;
using Artemis.Modifiers;

namespace Artemis.Simulation
{
    /// <summary>
    /// Cloth material properties
    /// </summary>
    public class ClothMaterial
    {
        public string Name { get; set; } = "Cloth";
        public float StretchStiffness { get; set; } = 1000f;      // Resistance to stretching
        public float BendStiffness { get; set; } = 50f;           // Resistance to bending
        public float ShearStiffness { get; set; } = 500f;         // Resistance to shearing
        public float Damping { get; set; } = 0.02f;               // Velocity damping
        public float Friction { get; set; } = 0.5f;               // Collision friction
        public float AirResistance { get; set; } = 0.02f;         // Air drag coefficient
        public float Density { get; set; } = 0.3f;                // kg/m² (mass per unit area)
        public float Thickness { get; set; } = 0.002f;            // meters
        public float TearThreshold { get; set; } = 10f;           // Force to tear (0 = unbreakable)
        public bool CanTear { get; set; } = true;

        // Presets
        public static ClothMaterial Cotton() => new()
        {
            Name = "Cotton",
            StretchStiffness = 800f,
            BendStiffness = 30f,
            ShearStiffness = 400f,
            Damping = 0.03f,
            AirResistance = 0.03f,
            Density = 0.2f,
            TearThreshold = 8f
        };

        public static ClothMaterial Silk() => new()
        {
            Name = "Silk",
            StretchStiffness = 600f,
            BendStiffness = 10f,
            ShearStiffness = 300f,
            Damping = 0.01f,
            AirResistance = 0.01f,
            Density = 0.1f,
            TearThreshold = 5f
        };

        public static ClothMaterial Denim() => new()
        {
            Name = "Denim",
            StretchStiffness = 2000f,
            BendStiffness = 100f,
            ShearStiffness = 1000f,
            Damping = 0.05f,
            AirResistance = 0.04f,
            Density = 0.5f,
            TearThreshold = 15f
        };

        public static ClothMaterial Leather() => new()
        {
            Name = "Leather",
            StretchStiffness = 5000f,
            BendStiffness = 200f,
            ShearStiffness = 2000f,
            Damping = 0.08f,
            AirResistance = 0.05f,
            Density = 0.8f,
            TearThreshold = 25f
        };

        public static ClothMaterial Rubber() => new()
        {
            Name = "Rubber",
            StretchStiffness = 300f,
            BendStiffness = 20f,
            ShearStiffness = 150f,
            Damping = 0.1f,
            AirResistance = 0.02f,
            Density = 1.2f,
            TearThreshold = 30f
        };

        public static ClothMaterial Flag() => new()
        {
            Name = "Flag",
            StretchStiffness = 1200f,
            BendStiffness = 5f,
            ShearStiffness = 600f,
            Damping = 0.01f,
            AirResistance = 0.05f,
            Density = 0.15f,
            TearThreshold = 20f,
            CanTear = false
        };

        public static ClothMaterial Sail() => new()
        {
            Name = "Sail",
            StretchStiffness = 3000f,
            BendStiffness = 80f,
            ShearStiffness = 1500f,
            Damping = 0.02f,
            AirResistance = 0.08f,
            Density = 0.4f,
            TearThreshold = 50f
        };

        public static ClothMaterial Curtain() => new()
        {
            Name = "Curtain",
            StretchStiffness = 500f,
            BendStiffness = 15f,
            ShearStiffness = 250f,
            Damping = 0.04f,
            AirResistance = 0.03f,
            Density = 0.25f,
            TearThreshold = 10f
        };
    }

    /// <summary>
    /// A single particle/vertex in the cloth mesh
    /// </summary>
    public class ClothParticle
    {
        public Vector3 Position { get; set; }
        public Vector3 PreviousPosition { get; set; }
        public Vector3 Velocity => Position - PreviousPosition;
        public Vector3 Acceleration { get; set; }
        public float Mass { get; set; } = 1f;
        public float InverseMass => IsPinned ? 0 : 1f / Mass;
        public bool IsPinned { get; set; }
        public bool IsActive { get; set; } = true;
        public Vector2 UV { get; set; }  // Texture coordinates
        public Vector3 Normal { get; set; } = Vector3.UnitY;
        public int GridX { get; set; }
        public int GridY { get; set; }
    }

    /// <summary>
    /// A spring constraint between two particles
    /// </summary>
    public class ClothConstraint
    {
        public ClothParticle ParticleA { get; set; }
        public ClothParticle ParticleB { get; set; }
        public float RestLength { get; set; }
        public float Stiffness { get; set; }
        public ConstraintType Type { get; set; }
        public bool IsBroken { get; set; }
        public float CurrentStretch => IsBroken ? 0 :
            ((ParticleA.Position - ParticleB.Position).Length() - RestLength) / RestLength;
    }

    public enum ConstraintType
    {
        Stretch,    // Horizontal/vertical neighbors
        Shear,      // Diagonal neighbors
        Bend        // Skip one neighbor (for bending resistance)
    }

    /// <summary>
    /// Result of a cloth collision
    /// </summary>
    public struct ClothCollision
    {
        public ClothParticle Particle;
        public Vector3 ContactPoint;
        public Vector3 ContactNormal;
        public float Penetration;
        public object Collider;
    }

    /// <summary>
    /// Spring-mass cloth simulation with force interaction
    /// </summary>
    public class ClothSimulation
    {
        private readonly List<ClothParticle> _particles = new();
        private readonly List<ClothConstraint> _constraints = new();
        private readonly List<Func<Vector3, Vector3>> _forceCallbacks = new();
        private readonly List<(Vector3 center, float radius, Func<Vector3, Vector3> getPosition)> _sphereColliders = new();
        private readonly List<(Vector3 min, Vector3 max)> _boxColliders = new();
        private readonly List<(Vector3 point, Vector3 normal)> _planeColliders = new();

        public ClothMaterial Material { get; set; } = new();
        public Vector3 Gravity { get; set; } = new(0, -9.81f, 0);
        public Vector3 Wind { get; set; } = Vector3.Zero;
        public float WindTurbulence { get; set; } = 0.3f;
        public int SolverIterations { get; set; } = 10;
        public int Width { get; private set; }
        public int Height { get; private set; }
        public float ParticleSpacing { get; private set; }
        public bool UseParallelUpdate { get; set; } = true;

        public IReadOnlyList<ClothParticle> Particles => _particles;
        public IReadOnlyList<ClothConstraint> Constraints => _constraints;
        public int ActiveParticleCount => _particles.Count(p => p.IsActive);
        public int BrokenConstraintCount => _constraints.Count(c => c.IsBroken);

        // Events
        public event Action<ClothConstraint>? OnConstraintBreak;
        public event Action<ClothCollision>? OnCollision;

        private Random _random = new();
        private float _time;

        /// <summary>
        /// Creates a rectangular cloth mesh
        /// </summary>
        /// <param name="origin">Top-left corner position</param>
        /// <param name="width">Number of particles horizontally</param>
        /// <param name="height">Number of particles vertically</param>
        /// <param name="spacing">Distance between particles</param>
        /// <param name="material">Cloth material properties</param>
        public ClothSimulation(Vector3 origin, int width, int height, float spacing, ClothMaterial? material = null)
        {
            Width = width;
            Height = height;
            ParticleSpacing = spacing;
            Material = material ?? new ClothMaterial();

            // Calculate mass per particle
            float areaPerParticle = spacing * spacing;
            float massPerParticle = Material.Density * areaPerParticle;

            // Create particles
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    var pos = origin + new Vector3(x * spacing, 0, y * spacing);
                    _particles.Add(new ClothParticle
                    {
                        Position = pos,
                        PreviousPosition = pos,
                        Mass = massPerParticle,
                        UV = new Vector2((float)x / (width - 1), (float)y / (height - 1)),
                        GridX = x,
                        GridY = y
                    });
                }
            }

            // Create constraints
            CreateConstraints();
        }

        private void CreateConstraints()
        {
            float diagLength = ParticleSpacing * MathF.Sqrt(2);

            for (int y = 0; y < Height; y++)
            {
                for (int x = 0; x < Width; x++)
                {
                    var p = GetParticle(x, y);
                    if (p == null) continue;

                    // Stretch constraints (horizontal and vertical)
                    if (x < Width - 1)
                    {
                        var right = GetParticle(x + 1, y);
                        if (right != null)
                            AddConstraint(p, right, ParticleSpacing, Material.StretchStiffness, ConstraintType.Stretch);
                    }

                    if (y < Height - 1)
                    {
                        var down = GetParticle(x, y + 1);
                        if (down != null)
                            AddConstraint(p, down, ParticleSpacing, Material.StretchStiffness, ConstraintType.Stretch);
                    }

                    // Shear constraints (diagonal)
                    if (x < Width - 1 && y < Height - 1)
                    {
                        var diagRight = GetParticle(x + 1, y + 1);
                        if (diagRight != null)
                            AddConstraint(p, diagRight, diagLength, Material.ShearStiffness, ConstraintType.Shear);
                    }

                    if (x > 0 && y < Height - 1)
                    {
                        var diagLeft = GetParticle(x - 1, y + 1);
                        if (diagLeft != null)
                            AddConstraint(p, diagLeft, diagLength, Material.ShearStiffness, ConstraintType.Shear);
                    }

                    // Bend constraints (skip one)
                    if (x < Width - 2)
                    {
                        var skip = GetParticle(x + 2, y);
                        if (skip != null)
                            AddConstraint(p, skip, ParticleSpacing * 2, Material.BendStiffness, ConstraintType.Bend);
                    }

                    if (y < Height - 2)
                    {
                        var skip = GetParticle(x, y + 2);
                        if (skip != null)
                            AddConstraint(p, skip, ParticleSpacing * 2, Material.BendStiffness, ConstraintType.Bend);
                    }
                }
            }
        }

        private void AddConstraint(ClothParticle a, ClothParticle b, float restLength, float stiffness, ConstraintType type)
        {
            _constraints.Add(new ClothConstraint
            {
                ParticleA = a,
                ParticleB = b,
                RestLength = restLength,
                Stiffness = stiffness,
                Type = type
            });
        }

        /// <summary>
        /// Get particle at grid position
        /// </summary>
        public ClothParticle? GetParticle(int x, int y)
        {
            if (x < 0 || x >= Width || y < 0 || y >= Height)
                return null;
            return _particles[y * Width + x];
        }

        /// <summary>
        /// Pin a particle (it won't move)
        /// </summary>
        public void PinParticle(int x, int y)
        {
            var p = GetParticle(x, y);
            if (p != null) p.IsPinned = true;
        }

        /// <summary>
        /// Unpin a particle
        /// </summary>
        public void UnpinParticle(int x, int y)
        {
            var p = GetParticle(x, y);
            if (p != null) p.IsPinned = false;
        }

        /// <summary>
        /// Pin the top edge
        /// </summary>
        public void PinTopEdge()
        {
            for (int x = 0; x < Width; x++)
                PinParticle(x, 0);
        }

        /// <summary>
        /// Pin corners only
        /// </summary>
        public void PinCorners()
        {
            PinParticle(0, 0);
            PinParticle(Width - 1, 0);
        }

        /// <summary>
        /// Move a pinned particle to a new position
        /// </summary>
        public void MovePinnedParticle(int x, int y, Vector3 newPosition)
        {
            var p = GetParticle(x, y);
            if (p != null && p.IsPinned)
            {
                p.PreviousPosition = p.Position;
                p.Position = newPosition;
            }
        }

        /// <summary>
        /// Add a custom force callback
        /// </summary>
        public void AddForce(Func<Vector3, Vector3> forceAtPosition)
        {
            _forceCallbacks.Add(forceAtPosition);
        }

        /// <summary>
        /// Apply force from a wind modifier
        /// </summary>
        public void ApplyWind(WindModifier wind)
        {
            Wind = new Vector3((float)wind.BaseDirection.X, (float)wind.BaseDirection.Y, (float)wind.BaseDirection.Z)
                   * (float)wind.BaseStrength;
            WindTurbulence = (float)wind.Turbulence;
        }

        /// <summary>
        /// Add sphere collider
        /// </summary>
        public void AddSphereCollider(Vector3 center, float radius, Func<Vector3, Vector3>? getPosition = null)
        {
            _sphereColliders.Add((center, radius, getPosition ?? (_ => center)));
        }

        /// <summary>
        /// Add box collider
        /// </summary>
        public void AddBoxCollider(Vector3 min, Vector3 max)
        {
            _boxColliders.Add((min, max));
        }

        /// <summary>
        /// Add plane collider
        /// </summary>
        public void AddPlaneCollider(Vector3 point, Vector3 normal)
        {
            _planeColliders.Add((point, Vector3.Normalize(normal)));
        }

        /// <summary>
        /// Update cloth simulation
        /// </summary>
        public void Update(float deltaTime)
        {
            _time += deltaTime;

            // Verlet integration with forces
            if (UseParallelUpdate && _particles.Count > 100)
            {
                Parallel.ForEach(_particles, particle =>
                {
                    if (!particle.IsActive || particle.IsPinned) return;
                    IntegrateParticle(particle, deltaTime);
                });
            }
            else
            {
                foreach (var particle in _particles)
                {
                    if (!particle.IsActive || particle.IsPinned) continue;
                    IntegrateParticle(particle, deltaTime);
                }
            }

            // Solve constraints
            for (int i = 0; i < SolverIterations; i++)
            {
                SolveConstraints();
            }

            // Handle collisions
            HandleCollisions();

            // Update normals
            UpdateNormals();
        }

        private void IntegrateParticle(ClothParticle p, float deltaTime)
        {
            // Calculate forces
            Vector3 force = Gravity * p.Mass;

            // Wind force (with turbulence)
            if (Wind.LengthSquared() > 0)
            {
                Vector3 windForce = Wind;
                windForce += new Vector3(
                    (float)(_random.NextDouble() - 0.5) * WindTurbulence,
                    (float)(_random.NextDouble() - 0.5) * WindTurbulence,
                    (float)(_random.NextDouble() - 0.5) * WindTurbulence
                ) * Wind.Length();

                // Wind affects cloth based on surface normal
                float windEffect = Math.Abs(Vector3.Dot(p.Normal, Vector3.Normalize(windForce)));
                force += windForce * windEffect * Material.AirResistance * 100;
            }

            // Air resistance (drag)
            Vector3 velocity = p.Position - p.PreviousPosition;
            force -= velocity * Material.AirResistance;

            // Custom forces
            foreach (var forceFunc in _forceCallbacks)
            {
                force += forceFunc(p.Position);
            }

            // Verlet integration
            p.Acceleration = force * p.InverseMass;

            Vector3 newPos = p.Position * 2 - p.PreviousPosition + p.Acceleration * deltaTime * deltaTime;
            newPos = Vector3.Lerp(newPos, p.Position + velocity, Material.Damping);

            p.PreviousPosition = p.Position;
            p.Position = newPos;
        }

        private void SolveConstraints()
        {
            foreach (var c in _constraints)
            {
                if (c.IsBroken) continue;

                var p1 = c.ParticleA;
                var p2 = c.ParticleB;

                if (!p1.IsActive || !p2.IsActive) continue;

                Vector3 delta = p2.Position - p1.Position;
                float currentLength = delta.Length();

                if (currentLength < 0.0001f) continue;

                float diff = (currentLength - c.RestLength) / currentLength;

                // Check for tearing
                if (Material.CanTear && c.Type == ConstraintType.Stretch)
                {
                    float strain = Math.Abs(currentLength - c.RestLength) / c.RestLength;
                    if (strain > Material.TearThreshold)
                    {
                        c.IsBroken = true;
                        OnConstraintBreak?.Invoke(c);
                        continue;
                    }
                }

                // Apply correction based on stiffness
                float stiffnessRatio = Math.Min(1f, c.Stiffness * 0.001f);
                Vector3 correction = delta * diff * 0.5f * stiffnessRatio;

                float w1 = p1.InverseMass;
                float w2 = p2.InverseMass;
                float totalW = w1 + w2;

                if (totalW > 0)
                {
                    p1.Position += correction * (w1 / totalW);
                    p2.Position -= correction * (w2 / totalW);
                }
            }
        }

        private void HandleCollisions()
        {
            foreach (var p in _particles)
            {
                if (!p.IsActive || p.IsPinned) continue;

                // Sphere colliders
                foreach (var (center, radius, getPos) in _sphereColliders)
                {
                    Vector3 spherePos = getPos(center);
                    Vector3 toParticle = p.Position - spherePos;
                    float dist = toParticle.Length();

                    if (dist < radius + Material.Thickness)
                    {
                        Vector3 normal = dist > 0.0001f ? toParticle / dist : Vector3.UnitY;
                        float penetration = radius + Material.Thickness - dist;

                        p.Position += normal * penetration;

                        // Friction
                        Vector3 velocity = p.Position - p.PreviousPosition;
                        Vector3 tangent = velocity - Vector3.Dot(velocity, normal) * normal;
                        p.PreviousPosition += tangent * Material.Friction;

                        OnCollision?.Invoke(new ClothCollision
                        {
                            Particle = p,
                            ContactPoint = spherePos + normal * radius,
                            ContactNormal = normal,
                            Penetration = penetration
                        });
                    }
                }

                // Plane colliders
                foreach (var (point, normal) in _planeColliders)
                {
                    float dist = Vector3.Dot(p.Position - point, normal);

                    if (dist < Material.Thickness)
                    {
                        float penetration = Material.Thickness - dist;
                        p.Position += normal * penetration;

                        // Friction
                        Vector3 velocity = p.Position - p.PreviousPosition;
                        Vector3 tangent = velocity - Vector3.Dot(velocity, normal) * normal;
                        p.PreviousPosition += tangent * Material.Friction;

                        OnCollision?.Invoke(new ClothCollision
                        {
                            Particle = p,
                            ContactPoint = p.Position - normal * Material.Thickness,
                            ContactNormal = normal,
                            Penetration = penetration
                        });
                    }
                }

                // Box colliders (simple AABB)
                foreach (var (min, max) in _boxColliders)
                {
                    if (p.Position.X > min.X && p.Position.X < max.X &&
                        p.Position.Y > min.Y && p.Position.Y < max.Y &&
                        p.Position.Z > min.Z && p.Position.Z < max.Z)
                    {
                        // Find closest face
                        float[] distances = {
                            p.Position.X - min.X, max.X - p.Position.X,
                            p.Position.Y - min.Y, max.Y - p.Position.Y,
                            p.Position.Z - min.Z, max.Z - p.Position.Z
                        };
                        Vector3[] normals = {
                            -Vector3.UnitX, Vector3.UnitX,
                            -Vector3.UnitY, Vector3.UnitY,
                            -Vector3.UnitZ, Vector3.UnitZ
                        };

                        int minIdx = 0;
                        for (int i = 1; i < 6; i++)
                        {
                            if (distances[i] < distances[minIdx])
                                minIdx = i;
                        }

                        p.Position += normals[minIdx] * (distances[minIdx] + Material.Thickness);

                        Vector3 velocity = p.Position - p.PreviousPosition;
                        Vector3 tangent = velocity - Vector3.Dot(velocity, normals[minIdx]) * normals[minIdx];
                        p.PreviousPosition += tangent * Material.Friction;
                    }
                }
            }
        }

        private void UpdateNormals()
        {
            // Reset normals
            foreach (var p in _particles)
            {
                p.Normal = Vector3.Zero;
            }

            // Calculate face normals and accumulate to vertices
            for (int y = 0; y < Height - 1; y++)
            {
                for (int x = 0; x < Width - 1; x++)
                {
                    var p00 = GetParticle(x, y);
                    var p10 = GetParticle(x + 1, y);
                    var p01 = GetParticle(x, y + 1);
                    var p11 = GetParticle(x + 1, y + 1);

                    if (p00 == null || p10 == null || p01 == null || p11 == null) continue;

                    // First triangle
                    Vector3 n1 = Vector3.Cross(p10.Position - p00.Position, p01.Position - p00.Position);
                    p00.Normal += n1;
                    p10.Normal += n1;
                    p01.Normal += n1;

                    // Second triangle
                    Vector3 n2 = Vector3.Cross(p01.Position - p11.Position, p10.Position - p11.Position);
                    p11.Normal += n2;
                    p01.Normal += n2;
                    p10.Normal += n2;
                }
            }

            // Normalize
            foreach (var p in _particles)
            {
                if (p.Normal.LengthSquared() > 0.0001f)
                    p.Normal = Vector3.Normalize(p.Normal);
                else
                    p.Normal = Vector3.UnitY;
            }
        }

        /// <summary>
        /// Apply an impulse at a point
        /// </summary>
        public void ApplyImpulse(Vector3 point, Vector3 impulse, float radius)
        {
            foreach (var p in _particles)
            {
                if (p.IsPinned) continue;

                float dist = (p.Position - point).Length();
                if (dist < radius)
                {
                    float factor = 1f - (dist / radius);
                    p.Position += impulse * factor * p.InverseMass;
                }
            }
        }

        /// <summary>
        /// Tear cloth at a line
        /// </summary>
        public void Tear(Vector3 start, Vector3 end, float width)
        {
            foreach (var c in _constraints)
            {
                if (c.IsBroken || c.Type != ConstraintType.Stretch) continue;

                Vector3 mid = (c.ParticleA.Position + c.ParticleB.Position) * 0.5f;

                // Distance from line
                Vector3 line = end - start;
                float t = Math.Clamp(Vector3.Dot(mid - start, line) / line.LengthSquared(), 0, 1);
                Vector3 closest = start + line * t;
                float dist = (mid - closest).Length();

                if (dist < width)
                {
                    c.IsBroken = true;
                    OnConstraintBreak?.Invoke(c);
                }
            }
        }

        /// <summary>
        /// Get vertices for rendering (triangle list)
        /// </summary>
        public List<(Vector3 pos, Vector3 normal, Vector2 uv)> GetRenderVertices()
        {
            var vertices = new List<(Vector3, Vector3, Vector2)>();

            for (int y = 0; y < Height - 1; y++)
            {
                for (int x = 0; x < Width - 1; x++)
                {
                    var p00 = GetParticle(x, y);
                    var p10 = GetParticle(x + 1, y);
                    var p01 = GetParticle(x, y + 1);
                    var p11 = GetParticle(x + 1, y + 1);

                    if (p00 == null || p10 == null || p01 == null || p11 == null) continue;

                    // Check if quad has broken constraints (skip torn areas)
                    bool hasBreak = _constraints.Any(c =>
                        c.IsBroken &&
                        ((c.ParticleA == p00 || c.ParticleA == p10 || c.ParticleA == p01 || c.ParticleA == p11) &&
                         (c.ParticleB == p00 || c.ParticleB == p10 || c.ParticleB == p01 || c.ParticleB == p11)));

                    if (hasBreak) continue;

                    // First triangle
                    vertices.Add((p00.Position, p00.Normal, p00.UV));
                    vertices.Add((p10.Position, p10.Normal, p10.UV));
                    vertices.Add((p01.Position, p01.Normal, p01.UV));

                    // Second triangle
                    vertices.Add((p10.Position, p10.Normal, p10.UV));
                    vertices.Add((p11.Position, p11.Normal, p11.UV));
                    vertices.Add((p01.Position, p01.Normal, p01.UV));
                }
            }

            return vertices;
        }

        /// <summary>
        /// Clear all colliders
        /// </summary>
        public void ClearColliders()
        {
            _sphereColliders.Clear();
            _boxColliders.Clear();
            _planeColliders.Clear();
        }

        /// <summary>
        /// Reset cloth to initial state
        /// </summary>
        public void Reset(Vector3 origin)
        {
            int i = 0;
            for (int y = 0; y < Height; y++)
            {
                for (int x = 0; x < Width; x++)
                {
                    var pos = origin + new Vector3(x * ParticleSpacing, 0, y * ParticleSpacing);
                    _particles[i].Position = pos;
                    _particles[i].PreviousPosition = pos;
                    _particles[i].Acceleration = Vector3.Zero;
                    i++;
                }
            }

            foreach (var c in _constraints)
            {
                c.IsBroken = false;
            }
        }
    }
}
