using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Collision;
using Artemis.Core;
using Artemis.Materials;
using Artemis.Particles;

namespace Artemis.Destruction
{
    /// <summary>
    /// Fracture pattern types for object destruction.
    /// </summary>
    public enum FracturePattern
    {
        /// <summary>Voronoi cell-based fracture (realistic).</summary>
        Voronoi,
        /// <summary>Radial fracture from impact point.</summary>
        Radial,
        /// <summary>Uniform grid-based fracture.</summary>
        Uniform,
        /// <summary>Random chunk fracture.</summary>
        Random,
        /// <summary>Shatter into many small pieces (glass-like).</summary>
        Shatter,
        /// <summary>Splinter pattern (wood-like).</summary>
        Splinter,
        /// <summary>Brick/mortar pattern.</summary>
        Brick
    }

    /// <summary>
    /// Configuration for fracture behavior.
    /// </summary>
    public class FractureConfig
    {
        /// <summary>Pattern to use for fracturing.</summary>
        public FracturePattern Pattern { get; set; } = FracturePattern.Voronoi;

        /// <summary>Minimum impulse required to cause fracture.</summary>
        public double FractureThreshold { get; set; } = 100.0;

        /// <summary>Minimum number of pieces to create.</summary>
        public int MinPieces { get; set; } = 3;

        /// <summary>Maximum number of pieces to create.</summary>
        public int MaxPieces { get; set; } = 20;

        /// <summary>Minimum fragment size as fraction of original.</summary>
        public double MinFragmentSize { get; set; } = 0.05;

        /// <summary>Whether to generate debris particles.</summary>
        public bool GenerateDebris { get; set; } = true;

        /// <summary>Number of debris particles per fracture.</summary>
        public int DebrisCount { get; set; } = 50;

        /// <summary>Whether fragments inherit velocity.</summary>
        public bool InheritVelocity { get; set; } = true;

        /// <summary>Additional velocity applied to fragments (explosion force).</summary>
        public double FragmentVelocityScale { get; set; } = 1.0;

        /// <summary>Whether fragments can further fracture.</summary>
        public bool AllowCascadeFracture { get; set; } = true;

        /// <summary>Depth limit for cascade fractures.</summary>
        public int MaxCascadeDepth { get; set; } = 3;
    }

    /// <summary>
    /// Result of a fracture operation.
    /// </summary>
    public class FractureResult
    {
        /// <summary>The original body that was fractured.</summary>
        public RigidBody OriginalBody { get; init; } = null!;

        /// <summary>List of fragment bodies created.</summary>
        public List<RigidBody> Fragments { get; } = new();

        /// <summary>Debris particles created.</summary>
        public List<Particle> Debris { get; } = new();

        /// <summary>Total mass of all fragments.</summary>
        public double TotalMass { get; set; }

        /// <summary>Point where fracture originated.</summary>
        public Vector3D ImpactPoint { get; init; }

        /// <summary>Direction of impact.</summary>
        public Vector3D ImpactDirection { get; init; }

        /// <summary>Cascade depth of this fracture.</summary>
        public int CascadeDepth { get; init; }
    }

    /// <summary>
    /// System for realistic object destruction and shattering.
    /// Supports various fracture patterns and debris generation.
    /// </summary>
    public class FractureSystem
    {
        #region Fields

        private readonly Random _random;
        private readonly Dictionary<RigidBody, FractureConfig> _fractureConfigs;
        private readonly List<FractureResult> _pendingFractures;
        private int _currentCascadeDepth;

        #endregion

        #region Properties

        /// <summary>Default fracture configuration.</summary>
        public FractureConfig DefaultConfig { get; set; } = new();

        /// <summary>Event raised when an object fractures.</summary>
        public event Action<FractureResult>? OnFracture;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new fracture system.
        /// </summary>
        public FractureSystem(int? seed = null)
        {
            _random = seed.HasValue ? new Random(seed.Value) : new Random();
            _fractureConfigs = new Dictionary<RigidBody, FractureConfig>();
            _pendingFractures = new List<FractureResult>();
        }

        #endregion

        #region Configuration

        /// <summary>
        /// Makes a body fracturable with the default config.
        /// </summary>
        public void MakeFracturable(RigidBody body)
        {
            MakeFracturable(body, DefaultConfig);
        }

        /// <summary>
        /// Makes a body fracturable with a custom config.
        /// </summary>
        public void MakeFracturable(RigidBody body, FractureConfig config)
        {
            _fractureConfigs[body] = config;
        }

        /// <summary>
        /// Gets whether a body is fracturable.
        /// </summary>
        public bool IsFracturable(RigidBody body)
        {
            return _fractureConfigs.ContainsKey(body);
        }

        /// <summary>
        /// Removes fracturability from a body.
        /// </summary>
        public void RemoveFracturable(RigidBody body)
        {
            _fractureConfigs.Remove(body);
        }

        #endregion

        #region Fracture Detection

        /// <summary>
        /// Checks if a collision should cause fracture and performs it.
        /// </summary>
        public FractureResult? CheckAndFracture(CollisionInfo collision, double impulse)
        {
            if (collision.BodyA is RigidBody rbA && IsFracturable(rbA))
            {
                var config = _fractureConfigs[rbA];
                if (impulse >= config.FractureThreshold)
                {
                    return Fracture(rbA, collision.ContactPoint, collision.Normal, impulse);
                }
            }

            if (collision.BodyB is RigidBody rbB && IsFracturable(rbB))
            {
                var config = _fractureConfigs[rbB];
                if (impulse >= config.FractureThreshold)
                {
                    return Fracture(rbB, collision.ContactPoint, -collision.Normal, impulse);
                }
            }

            return null;
        }

        /// <summary>
        /// Fractures a body at a specific point with given impact direction.
        /// </summary>
        public FractureResult Fracture(RigidBody body, Vector3D impactPoint, Vector3D impactDirection, double force)
        {
            var config = _fractureConfigs.GetValueOrDefault(body) ?? DefaultConfig;

            var result = new FractureResult
            {
                OriginalBody = body,
                ImpactPoint = impactPoint,
                ImpactDirection = impactDirection.Normalized,
                CascadeDepth = _currentCascadeDepth
            };

            // Generate fragments based on pattern
            switch (config.Pattern)
            {
                case FracturePattern.Voronoi:
                    GenerateVoronoiFragments(body, impactPoint, config, result);
                    break;
                case FracturePattern.Radial:
                    GenerateRadialFragments(body, impactPoint, config, result);
                    break;
                case FracturePattern.Uniform:
                    GenerateUniformFragments(body, config, result);
                    break;
                case FracturePattern.Random:
                    GenerateRandomFragments(body, config, result);
                    break;
                case FracturePattern.Shatter:
                    GenerateShatterFragments(body, impactPoint, config, result);
                    break;
                case FracturePattern.Splinter:
                    GenerateSplinterFragments(body, impactPoint, impactDirection, config, result);
                    break;
                case FracturePattern.Brick:
                    GenerateBrickFragments(body, config, result);
                    break;
            }

            // Apply velocities to fragments
            if (config.InheritVelocity)
            {
                foreach (var fragment in result.Fragments)
                {
                    // Inherit base velocity
                    fragment.Velocity = body.Velocity;
                    fragment.AngularVelocity = body.AngularVelocity;

                    // Add explosion velocity from impact point
                    var toFragment = (fragment.Position - impactPoint).Normalized;
                    var explosionVel = toFragment * force * config.FragmentVelocityScale * 0.01;
                    fragment.Velocity += explosionVel;

                    // Add random angular velocity
                    fragment.AngularVelocity += new Vector3D(
                        (_random.NextDouble() - 0.5) * 5,
                        (_random.NextDouble() - 0.5) * 5,
                        (_random.NextDouble() - 0.5) * 5
                    );
                }
            }

            // Generate debris particles
            if (config.GenerateDebris)
            {
                GenerateDebris(body, impactPoint, impactDirection, force, config, result);
            }

            // Calculate total mass
            result.TotalMass = 0;
            foreach (var fragment in result.Fragments)
            {
                result.TotalMass += fragment.Mass;
            }

            // Remove original from fracturable tracking
            _fractureConfigs.Remove(body);

            // Add fragments as fracturable if cascade allowed
            if (config.AllowCascadeFracture && _currentCascadeDepth < config.MaxCascadeDepth)
            {
                var cascadeConfig = new FractureConfig
                {
                    Pattern = config.Pattern,
                    FractureThreshold = config.FractureThreshold * 0.7, // Lower threshold for fragments
                    MinPieces = Math.Max(2, config.MinPieces / 2),
                    MaxPieces = Math.Max(3, config.MaxPieces / 2),
                    MinFragmentSize = config.MinFragmentSize,
                    GenerateDebris = config.GenerateDebris,
                    DebrisCount = config.DebrisCount / 2,
                    InheritVelocity = config.InheritVelocity,
                    FragmentVelocityScale = config.FragmentVelocityScale * 0.5,
                    AllowCascadeFracture = true,
                    MaxCascadeDepth = config.MaxCascadeDepth
                };

                foreach (var fragment in result.Fragments)
                {
                    if (fragment.Mass > body.Mass * config.MinFragmentSize * 2)
                    {
                        _fractureConfigs[fragment] = cascadeConfig;
                    }
                }
            }

            OnFracture?.Invoke(result);
            return result;
        }

        #endregion

        #region Fragment Generation

        private void GenerateVoronoiFragments(RigidBody body, Vector3D impactPoint, FractureConfig config, FractureResult result)
        {
            int pieceCount = _random.Next(config.MinPieces, config.MaxPieces + 1);
            var halfExtents = body.HalfExtents;
            double totalVolume = halfExtents.X * halfExtents.Y * halfExtents.Z * 8;

            // Generate Voronoi seed points
            var seeds = new List<Vector3D>();
            seeds.Add(body.Transform.InverseTransformPoint(impactPoint)); // Impact point is always a seed

            for (int i = 1; i < pieceCount; i++)
            {
                var localSeed = new Vector3D(
                    (_random.NextDouble() * 2 - 1) * halfExtents.X,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Y,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Z
                );
                seeds.Add(localSeed);
            }

            // Create fragments (simplified - actual Voronoi would be more complex)
            double massPerFragment = body.Mass / pieceCount;
            double sizeScale = Math.Pow(1.0 / pieceCount, 1.0 / 3.0);

            for (int i = 0; i < pieceCount; i++)
            {
                var worldPos = body.Transform.TransformPoint(seeds[i]);
                var fragmentSize = halfExtents * sizeScale * (0.5 + _random.NextDouble() * 0.5);

                // Ensure minimum size
                fragmentSize = new Vector3D(
                    Math.Max(fragmentSize.X, halfExtents.X * config.MinFragmentSize),
                    Math.Max(fragmentSize.Y, halfExtents.Y * config.MinFragmentSize),
                    Math.Max(fragmentSize.Z, halfExtents.Z * config.MinFragmentSize)
                );

                var fragment = RigidBody.CreateBox(worldPos, fragmentSize, massPerFragment, body.Material);
                fragment.Rotation = body.Rotation;
                result.Fragments.Add(fragment);
            }
        }

        private void GenerateRadialFragments(RigidBody body, Vector3D impactPoint, FractureConfig config, FractureResult result)
        {
            int pieceCount = _random.Next(config.MinPieces, config.MaxPieces + 1);
            var halfExtents = body.HalfExtents;
            double massPerFragment = body.Mass / pieceCount;

            // Create radial slices around impact point
            double angleStep = 2 * Math.PI / pieceCount;
            double sizeScale = Math.Pow(1.0 / pieceCount, 1.0 / 3.0);

            var localImpact = body.Transform.InverseTransformPoint(impactPoint);

            for (int i = 0; i < pieceCount; i++)
            {
                double angle = i * angleStep;
                var offset = new Vector3D(
                    Math.Cos(angle) * halfExtents.X * 0.5,
                    0,
                    Math.Sin(angle) * halfExtents.Z * 0.5
                );

                var fragmentPos = body.Position + body.Transform.TransformDirection(offset);
                var fragmentSize = halfExtents * sizeScale;

                var fragment = RigidBody.CreateBox(fragmentPos, fragmentSize, massPerFragment, body.Material);
                fragment.Rotation = body.Rotation;
                result.Fragments.Add(fragment);
            }
        }

        private void GenerateUniformFragments(RigidBody body, FractureConfig config, FractureResult result)
        {
            // Calculate grid dimensions
            int totalPieces = _random.Next(config.MinPieces, config.MaxPieces + 1);
            int gridSize = Math.Max(1, (int)Math.Ceiling(Math.Pow(totalPieces, 1.0 / 3.0)));

            var halfExtents = body.HalfExtents;
            var fragmentSize = new Vector3D(
                halfExtents.X / gridSize,
                halfExtents.Y / gridSize,
                halfExtents.Z / gridSize
            );

            double massPerFragment = body.Mass / (gridSize * gridSize * gridSize);

            for (int x = 0; x < gridSize; x++)
            {
                for (int y = 0; y < gridSize; y++)
                {
                    for (int z = 0; z < gridSize; z++)
                    {
                        var localOffset = new Vector3D(
                            -halfExtents.X + fragmentSize.X * (2 * x + 1),
                            -halfExtents.Y + fragmentSize.Y * (2 * y + 1),
                            -halfExtents.Z + fragmentSize.Z * (2 * z + 1)
                        );

                        var fragmentPos = body.Position + body.Transform.TransformDirection(localOffset);
                        var fragment = RigidBody.CreateBox(fragmentPos, fragmentSize, massPerFragment, body.Material);
                        fragment.Rotation = body.Rotation;
                        result.Fragments.Add(fragment);
                    }
                }
            }
        }

        private void GenerateRandomFragments(RigidBody body, FractureConfig config, FractureResult result)
        {
            int pieceCount = _random.Next(config.MinPieces, config.MaxPieces + 1);
            var halfExtents = body.HalfExtents;
            double massPerFragment = body.Mass / pieceCount;

            for (int i = 0; i < pieceCount; i++)
            {
                var localOffset = new Vector3D(
                    (_random.NextDouble() * 2 - 1) * halfExtents.X * 0.8,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Y * 0.8,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Z * 0.8
                );

                var fragmentSize = new Vector3D(
                    halfExtents.X * (config.MinFragmentSize + _random.NextDouble() * (0.5 - config.MinFragmentSize)),
                    halfExtents.Y * (config.MinFragmentSize + _random.NextDouble() * (0.5 - config.MinFragmentSize)),
                    halfExtents.Z * (config.MinFragmentSize + _random.NextDouble() * (0.5 - config.MinFragmentSize))
                );

                var fragmentPos = body.Position + body.Transform.TransformDirection(localOffset);
                var fragment = RigidBody.CreateBox(fragmentPos, fragmentSize, massPerFragment, body.Material);
                fragment.Rotation = body.Rotation;
                result.Fragments.Add(fragment);
            }
        }

        private void GenerateShatterFragments(RigidBody body, Vector3D impactPoint, FractureConfig config, FractureResult result)
        {
            // Glass-like shatter - many small pieces near impact, larger pieces away
            int pieceCount = config.MaxPieces * 2; // More pieces for shatter
            var halfExtents = body.HalfExtents;
            var localImpact = body.Transform.InverseTransformPoint(impactPoint);
            double maxDist = halfExtents.Magnitude;

            for (int i = 0; i < pieceCount; i++)
            {
                var localOffset = new Vector3D(
                    (_random.NextDouble() * 2 - 1) * halfExtents.X,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Y,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Z
                );

                double distFromImpact = (localOffset - localImpact).Magnitude / maxDist;

                // Size based on distance from impact (smaller near impact)
                double sizeMultiplier = 0.1 + distFromImpact * 0.4;
                var fragmentSize = halfExtents * sizeMultiplier * config.MinFragmentSize * 3;

                double massPerFragment = body.Mass / pieceCount;

                var fragmentPos = body.Position + body.Transform.TransformDirection(localOffset);
                var fragment = RigidBody.CreateBox(fragmentPos, fragmentSize, massPerFragment, body.Material);
                fragment.Rotation = body.Rotation;
                result.Fragments.Add(fragment);
            }
        }

        private void GenerateSplinterFragments(RigidBody body, Vector3D impactPoint, Vector3D impactDir, FractureConfig config, FractureResult result)
        {
            // Wood-like splintering - elongated pieces along grain direction
            int pieceCount = _random.Next(config.MinPieces, config.MaxPieces + 1);
            var halfExtents = body.HalfExtents;
            double massPerFragment = body.Mass / pieceCount;

            // Determine grain direction (longest axis)
            Vector3D grainDir;
            if (halfExtents.Y > halfExtents.X && halfExtents.Y > halfExtents.Z)
                grainDir = Vector3D.Up;
            else if (halfExtents.X > halfExtents.Z)
                grainDir = Vector3D.Right;
            else
                grainDir = Vector3D.Forward;

            grainDir = body.Transform.TransformDirection(grainDir);

            for (int i = 0; i < pieceCount; i++)
            {
                var localOffset = new Vector3D(
                    (_random.NextDouble() * 2 - 1) * halfExtents.X,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Y,
                    (_random.NextDouble() * 2 - 1) * halfExtents.Z
                );

                // Elongated along grain
                var fragmentSize = new Vector3D(
                    halfExtents.X * config.MinFragmentSize * (1 + _random.NextDouble()),
                    halfExtents.Y * (0.3 + _random.NextDouble() * 0.5),
                    halfExtents.Z * config.MinFragmentSize * (1 + _random.NextDouble())
                );

                var fragmentPos = body.Position + body.Transform.TransformDirection(localOffset);
                var fragment = RigidBody.CreateBox(fragmentPos, fragmentSize, massPerFragment, body.Material);
                fragment.Rotation = body.Rotation;
                result.Fragments.Add(fragment);
            }
        }

        private void GenerateBrickFragments(RigidBody body, FractureConfig config, FractureResult result)
        {
            // Brick/mortar pattern - regular rectangular pieces
            int gridX = _random.Next(2, 5);
            int gridY = _random.Next(2, 4);
            int gridZ = _random.Next(2, 5);

            var halfExtents = body.HalfExtents;
            var fragmentSize = new Vector3D(
                halfExtents.X / gridX,
                halfExtents.Y / gridY,
                halfExtents.Z / gridZ
            );

            int totalPieces = gridX * gridY * gridZ;
            double massPerFragment = body.Mass / totalPieces;

            for (int x = 0; x < gridX; x++)
            {
                for (int y = 0; y < gridY; y++)
                {
                    for (int z = 0; z < gridZ; z++)
                    {
                        // Add small random offset for natural look
                        var jitter = new Vector3D(
                            (_random.NextDouble() - 0.5) * fragmentSize.X * 0.1,
                            (_random.NextDouble() - 0.5) * fragmentSize.Y * 0.1,
                            (_random.NextDouble() - 0.5) * fragmentSize.Z * 0.1
                        );

                        var localOffset = new Vector3D(
                            -halfExtents.X + fragmentSize.X * (2 * x + 1),
                            -halfExtents.Y + fragmentSize.Y * (2 * y + 1),
                            -halfExtents.Z + fragmentSize.Z * (2 * z + 1)
                        ) + jitter;

                        var fragmentPos = body.Position + body.Transform.TransformDirection(localOffset);
                        var fragment = RigidBody.CreateBox(fragmentPos, fragmentSize * 0.95, massPerFragment, body.Material);
                        fragment.Rotation = body.Rotation;
                        result.Fragments.Add(fragment);
                    }
                }
            }
        }

        private void GenerateDebris(RigidBody body, Vector3D impactPoint, Vector3D impactDir, double force, FractureConfig config, FractureResult result)
        {
            for (int i = 0; i < config.DebrisCount; i++)
            {
                // Random direction biased toward impact direction
                var dir = new Vector3D(
                    _random.NextDouble() * 2 - 1,
                    _random.NextDouble() * 2 - 1,
                    _random.NextDouble() * 2 - 1
                ).Normalized;
                dir = (dir + impactDir * 0.5).Normalized;

                var velocity = dir * force * 0.02 * (0.5 + _random.NextDouble());

                var debris = new Particle
                {
                    Position = impactPoint + dir * 0.1,
                    Velocity = body.Velocity + velocity,
                    Size = 0.02 + _random.NextDouble() * 0.05,
                    Mass = 0.01,
                    Lifetime = 2.0 + _random.NextDouble() * 3.0,
                    Color = body.Material?.Color ?? 0xFF808080
                };

                result.Debris.Add(debris);
            }
        }

        #endregion

        #region Preset Configurations

        /// <summary>Creates a glass fracture configuration.</summary>
        public static FractureConfig Glass() => new()
        {
            Pattern = FracturePattern.Shatter,
            FractureThreshold = 30,
            MinPieces = 10,
            MaxPieces = 50,
            MinFragmentSize = 0.02,
            GenerateDebris = true,
            DebrisCount = 100,
            FragmentVelocityScale = 2.0,
            AllowCascadeFracture = false
        };

        /// <summary>Creates a wood fracture configuration.</summary>
        public static FractureConfig Wood() => new()
        {
            Pattern = FracturePattern.Splinter,
            FractureThreshold = 150,
            MinPieces = 3,
            MaxPieces = 10,
            MinFragmentSize = 0.1,
            GenerateDebris = true,
            DebrisCount = 30,
            FragmentVelocityScale = 0.5,
            AllowCascadeFracture = true,
            MaxCascadeDepth = 2
        };

        /// <summary>Creates a stone/concrete fracture configuration.</summary>
        public static FractureConfig Stone() => new()
        {
            Pattern = FracturePattern.Voronoi,
            FractureThreshold = 300,
            MinPieces = 4,
            MaxPieces = 12,
            MinFragmentSize = 0.15,
            GenerateDebris = true,
            DebrisCount = 50,
            FragmentVelocityScale = 0.3,
            AllowCascadeFracture = true,
            MaxCascadeDepth = 2
        };

        /// <summary>Creates a brick fracture configuration.</summary>
        public static FractureConfig Brick() => new()
        {
            Pattern = FracturePattern.Brick,
            FractureThreshold = 200,
            MinPieces = 6,
            MaxPieces = 20,
            MinFragmentSize = 0.1,
            GenerateDebris = true,
            DebrisCount = 40,
            FragmentVelocityScale = 0.4,
            AllowCascadeFracture = false
        };

        /// <summary>Creates a metal fracture configuration.</summary>
        public static FractureConfig Metal() => new()
        {
            Pattern = FracturePattern.Radial,
            FractureThreshold = 500,
            MinPieces = 2,
            MaxPieces = 6,
            MinFragmentSize = 0.2,
            GenerateDebris = true,
            DebrisCount = 20,
            FragmentVelocityScale = 0.2,
            AllowCascadeFracture = false
        };

        /// <summary>Creates an ice fracture configuration.</summary>
        public static FractureConfig Ice() => new()
        {
            Pattern = FracturePattern.Shatter,
            FractureThreshold = 20,
            MinPieces = 15,
            MaxPieces = 40,
            MinFragmentSize = 0.03,
            GenerateDebris = true,
            DebrisCount = 80,
            FragmentVelocityScale = 1.5,
            AllowCascadeFracture = false
        };

        #endregion
    }

    /// <summary>
    /// Represents an object made of particles that can erode/dissolve.
    /// Useful for sand sculptures, snow, dirt, etc.
    /// </summary>
    public class ErodibleBody
    {
        #region Fields

#if NETSTANDARD2_1
        private static readonly Random SharedRandom = new Random();
        private static readonly object SharedRandomLock = new object();
#endif

        private readonly List<Particle> _particles;
        private readonly SpatialHash _spatialHash;
        private readonly double _cohesion;
        private readonly double _erosionRate;

        #endregion

        #region Properties

        /// <summary>Gets all particles in this body.</summary>
        public IReadOnlyList<Particle> Particles => _particles;

        /// <summary>Gets the center of mass.</summary>
        public Vector3D CenterOfMass
        {
            get
            {
                if (_particles.Count == 0) return Vector3D.Zero;
                var sum = Vector3D.Zero;
                foreach (var p in _particles)
                    sum += p.Position;
                return sum / _particles.Count;
            }
        }

        /// <summary>Gets the total mass.</summary>
        public double TotalMass
        {
            get
            {
                double mass = 0;
                foreach (var p in _particles)
                    mass += p.Mass;
                return mass;
            }
        }

        #endregion

        #region Constructors

        /// <summary>
        /// Creates an erodible body from a list of particles.
        /// </summary>
        public ErodibleBody(IEnumerable<Particle> particles, double cohesion = 0.5, double erosionRate = 0.1)
        {
            _particles = new List<Particle>(particles);
            _cohesion = cohesion;
            _erosionRate = erosionRate;
            _spatialHash = new SpatialHash(0.2);
        }

        /// <summary>
        /// Creates an erodible body filling a box shape.
        /// </summary>
        public static ErodibleBody CreateBox(Vector3D center, Vector3D halfExtents, double particleSize, double cohesion = 0.5)
        {
            var particles = new List<Particle>();
            double spacing = particleSize * 2;

            int countX = (int)(halfExtents.X * 2 / spacing);
            int countY = (int)(halfExtents.Y * 2 / spacing);
            int countZ = (int)(halfExtents.Z * 2 / spacing);

            for (int x = 0; x < countX; x++)
            {
                for (int y = 0; y < countY; y++)
                {
                    for (int z = 0; z < countZ; z++)
                    {
                        var pos = center + new Vector3D(
                            -halfExtents.X + spacing * (x + 0.5),
                            -halfExtents.Y + spacing * (y + 0.5),
                            -halfExtents.Z + spacing * (z + 0.5)
                        );

                        particles.Add(new Particle
                        {
                            Position = pos,
                            Size = particleSize,
                            Mass = particleSize * particleSize * particleSize * 2000,
                            Lifetime = double.MaxValue,
                            Color = 0xFFC2B280 // Sand color
                        });
                    }
                }
            }

            return new ErodibleBody(particles, cohesion);
        }

        /// <summary>
        /// Creates an erodible sphere.
        /// </summary>
        public static ErodibleBody CreateSphere(Vector3D center, double radius, double particleSize, double cohesion = 0.5)
        {
            var particles = new List<Particle>();
            double spacing = particleSize * 2;
            double radiusSq = radius * radius;

            int count = (int)(radius * 2 / spacing) + 1;

            for (int x = 0; x < count; x++)
            {
                for (int y = 0; y < count; y++)
                {
                    for (int z = 0; z < count; z++)
                    {
                        var offset = new Vector3D(
                            -radius + spacing * x,
                            -radius + spacing * y,
                            -radius + spacing * z
                        );

                        if (offset.MagnitudeSquared <= radiusSq)
                        {
                            particles.Add(new Particle
                            {
                                Position = center + offset,
                                Size = particleSize,
                                Mass = particleSize * particleSize * particleSize * 2000,
                                Lifetime = double.MaxValue,
                                Color = 0xFFC2B280
                            });
                        }
                    }
                }
            }

            return new ErodibleBody(particles, cohesion);
        }

        #endregion

        #region Erosion

        private static double NextSharedDouble()
        {
#if NETSTANDARD2_1
            lock (SharedRandomLock)
            {
                return SharedRandom.NextDouble();
            }
#else
            return Random.Shared.NextDouble();
#endif
        }

        /// <summary>
        /// Applies wind erosion to the body.
        /// </summary>
        public List<Particle> ApplyWind(Vector3D windDirection, double windStrength, double deltaTime)
        {
            var erodedParticles = new List<Particle>();
            windDirection = windDirection.Normalized;

            for (int i = _particles.Count - 1; i >= 0; i--)
            {
                var p = _particles[i];

                // Check if particle is exposed (on surface)
                bool isExposed = IsParticleExposed(p, windDirection);

                if (isExposed)
                {
                    // Calculate erosion chance based on wind strength and cohesion
                    double erosionChance = windStrength * _erosionRate * (1 - _cohesion) * deltaTime;

                    if (NextSharedDouble() < erosionChance)
                    {
                        // Remove particle and add to eroded list
                        _particles.RemoveAt(i);

                        // Give it velocity in wind direction
                        p.Velocity = windDirection * windStrength * (0.5 + NextSharedDouble());
                        p.Lifetime = 5 + NextSharedDouble() * 5;
                        erodedParticles.Add(p);
                    }
                }
            }

            return erodedParticles;
        }

        /// <summary>
        /// Applies impact erosion (e.g., projectile hit).
        /// </summary>
        public List<Particle> ApplyImpact(Vector3D impactPoint, double impactForce, double impactRadius)
        {
            var erodedParticles = new List<Particle>();
            double radiusSq = impactRadius * impactRadius;

            for (int i = _particles.Count - 1; i >= 0; i--)
            {
                var p = _particles[i];
                var toParticle = p.Position - impactPoint;
                double distSq = toParticle.MagnitudeSquared;

                if (distSq < radiusSq)
                {
                    double dist = Math.Sqrt(distSq);
                    double forceFactor = 1 - (dist / impactRadius);
                    double erosionChance = forceFactor * impactForce * (1 - _cohesion) * 0.01;

                    if (NextSharedDouble() < erosionChance)
                    {
                        _particles.RemoveAt(i);

                        var direction = toParticle.Normalized;
                        p.Velocity = direction * impactForce * forceFactor * 0.1;
                        p.Lifetime = 3 + NextSharedDouble() * 3;
                        erodedParticles.Add(p);
                    }
                }
            }

            return erodedParticles;
        }

        private bool IsParticleExposed(Particle p, Vector3D checkDirection)
        {
            // Check if there's a neighbor particle in the check direction
            double checkDist = p.Size * 2.5;
            var checkPos = p.Position + checkDirection * checkDist;

            foreach (var other in _particles)
            {
                if (other == p) continue;
                if ((other.Position - checkPos).MagnitudeSquared < p.Size * p.Size * 4)
                {
                    return false;
                }
            }
            return true;
        }

        #endregion
    }
}
