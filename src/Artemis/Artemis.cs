// Artemis Physics Engine
// A lightweight, extensible physics engine for C# applications
// Compatible with Unity, MonoGame, and other frameworks

// Core types - re-exported for convenience
global using Artemis.Core;
global using Artemis.Bodies;
global using Artemis.Forces;
global using Artemis.Materials;
global using Artemis.Particles;
global using Artemis.Collision;
global using Artemis.Simulation;
global using Artemis.Compute;
global using Artemis.Destruction;
global using Artemis.Modifiers;

namespace Artemis
{
    /// <summary>
    /// Artemis Physics Engine - Main entry point and factory methods.
    /// </summary>
    public static class Physics
    {
        /// <summary>
        /// Gets the version of the Artemis Physics Engine.
        /// </summary>
        public static string Version => "1.0.0";

        /// <summary>
        /// Creates a new physics world with default settings.
        /// </summary>
        /// <returns>A new PhysicsWorld instance.</returns>
        public static PhysicsWorld CreateWorld()
        {
            return new PhysicsWorld();
        }

        /// <summary>
        /// Creates a new physics world with custom gravity.
        /// </summary>
        /// <param name="gravity">The gravity vector.</param>
        /// <returns>A new PhysicsWorld instance.</returns>
        public static PhysicsWorld CreateWorld(Vector3D gravity)
        {
            return new PhysicsWorld(gravity);
        }

        /// <summary>
        /// Creates a new physics world configured for Earth gravity.
        /// </summary>
        public static PhysicsWorld CreateEarthWorld()
        {
            return new PhysicsWorld(new Vector3D(0, -PhysicsConstants.EarthGravity, 0));
        }

        /// <summary>
        /// Creates a new physics world configured for Moon gravity.
        /// </summary>
        public static PhysicsWorld CreateMoonWorld()
        {
            return new PhysicsWorld(new Vector3D(0, -PhysicsConstants.MoonGravity, 0));
        }

        /// <summary>
        /// Creates a new physics world configured for Mars gravity.
        /// </summary>
        public static PhysicsWorld CreateMarsWorld()
        {
            return new PhysicsWorld(new Vector3D(0, -PhysicsConstants.MarsGravity, 0));
        }

        /// <summary>
        /// Creates a new physics world with zero gravity.
        /// </summary>
        public static PhysicsWorld CreateZeroGWorld()
        {
            return new PhysicsWorld(Vector3D.Zero);
        }

        /// <summary>
        /// Creates a new particle system.
        /// </summary>
        /// <param name="maxParticles">Maximum number of particles.</param>
        /// <returns>A new ParticleSystem instance.</returns>
        public static ParticleSystem CreateParticleSystem(int maxParticles = 10000)
        {
            return new ParticleSystem(maxParticles);
        }

        /// <summary>
        /// Creates a new sand simulation.
        /// </summary>
        /// <param name="bounds">Simulation bounds.</param>
        /// <param name="cellSize">Spatial grid cell size.</param>
        /// <returns>A new SandSimulation instance.</returns>
        public static SandSimulation CreateSandSimulation(AABB bounds, double cellSize = 0.5)
        {
            return new SandSimulation(bounds, cellSize);
        }

        #region Quick Create Methods

        /// <summary>
        /// Creates a dynamic sphere body.
        /// </summary>
        public static RigidBody CreateSphere(
            Vector3D position,
            double radius,
            double mass,
            PhysicsMaterial? material = null)
        {
            return RigidBody.CreateSphere(position, radius, mass, material);
        }

        /// <summary>
        /// Creates a dynamic box body.
        /// </summary>
        public static RigidBody CreateBox(
            Vector3D position,
            Vector3D halfExtents,
            double mass,
            PhysicsMaterial? material = null)
        {
            return RigidBody.CreateBox(position, halfExtents, mass, material);
        }

        /// <summary>
        /// Creates a static box body (for floors, walls, etc.).
        /// </summary>
        public static RigidBody CreateStaticBox(
            Vector3D position,
            Vector3D halfExtents,
            PhysicsMaterial? material = null)
        {
            return RigidBody.CreateStaticBox(position, halfExtents, material);
        }

        /// <summary>
        /// Creates a static sphere body.
        /// </summary>
        public static RigidBody CreateStaticSphere(
            Vector3D position,
            double radius,
            PhysicsMaterial? material = null)
        {
            return RigidBody.CreateStaticSphere(position, radius, material);
        }

        #endregion

        #region Force Factory Methods

        /// <summary>
        /// Creates a uniform gravity force.
        /// </summary>
        public static GravityForce CreateGravity(double magnitude = PhysicsConstants.EarthGravity)
        {
            return new GravityForce(magnitude);
        }

        /// <summary>
        /// Creates a point gravity (attraction to a point).
        /// </summary>
        public static PointGravityForce CreatePointGravity(
            Vector3D position,
            double mass,
            double? gravitationalConstant = null)
        {
            return gravitationalConstant.HasValue
                ? new PointGravityForce(position, mass, gravitationalConstant.Value)
                : new PointGravityForce(position, mass);
        }

        /// <summary>
        /// Creates an aerodynamic drag force.
        /// </summary>
        public static DragForce CreateDrag(
            double dragCoefficient = 0.47,
            double crossSectionArea = 1.0,
            double fluidDensity = 1.225)
        {
            return new DragForce(dragCoefficient, crossSectionArea, fluidDensity);
        }

        /// <summary>
        /// Creates a spring force.
        /// </summary>
        public static SpringForce CreateSpring(
            Vector3D anchor,
            double stiffness = 100.0,
            double restLength = 1.0,
            double damping = 0.1)
        {
            return new SpringForce(anchor, stiffness, restLength, damping);
        }

        /// <summary>
        /// Creates a repulsion force.
        /// </summary>
        public static RepulsionForce CreateRepulsion(
            Vector3D position,
            double strength = 100.0,
            double maxRange = 10.0)
        {
            return new RepulsionForce(position, strength, maxRange);
        }

        /// <summary>
        /// Creates a wind force.
        /// </summary>
        public static WindForce CreateWind(
            Vector3D direction,
            double strength = 10.0,
            double turbulence = 0.2)
        {
            return new WindForce(direction, strength, turbulence);
        }

        /// <summary>
        /// Creates a magnetic force.
        /// </summary>
        public static MagneticForce CreateMagnet(
            Vector3D position,
            double strength = 100.0,
            bool attracting = true)
        {
            return new MagneticForce(position, strength).Attracting(attracting);
        }

        /// <summary>
        /// Creates a buoyancy force for water simulation.
        /// </summary>
        public static BuoyancyForce CreateBuoyancy(
            double surfaceHeight = 0,
            double fluidDensity = 1000.0)
        {
            return new BuoyancyForce(surfaceHeight, fluidDensity);
        }

        /// <summary>
        /// Creates a vortex/whirlpool force.
        /// </summary>
        public static VortexForce CreateVortex(
            Vector3D center,
            double rotationStrength = 50.0,
            double pullStrength = 20.0)
        {
            return new VortexForce(center, null, rotationStrength, pullStrength);
        }

        /// <summary>
        /// Creates an explosion force.
        /// </summary>
        public static ExplosionForce CreateExplosion(
            Vector3D center,
            double force = 1000.0,
            double radius = 10.0)
        {
            return new ExplosionForce(center, force, radius);
        }

        /// <summary>
        /// Creates a friction force for a surface.
        /// </summary>
        public static FrictionForce CreateFriction(
            double staticCoefficient = 0.5,
            double kineticCoefficient = 0.3,
            double surfaceHeight = 0)
        {
            return new FrictionForce(staticCoefficient, kineticCoefficient, surfaceHeight);
        }

        #endregion

        #region Fluid Simulation

        /// <summary>
        /// Creates a fluid (SPH) simulation.
        /// </summary>
        public static FluidSimulation CreateFluidSimulation(AABB bounds, double smoothingRadius = 0.2)
        {
            return new FluidSimulation(bounds, smoothingRadius);
        }

        /// <summary>
        /// Creates a water simulation.
        /// </summary>
        public static FluidSimulation CreateWater(AABB bounds)
        {
            return FluidSimulation.Water(bounds);
        }

        /// <summary>
        /// Creates an oil simulation.
        /// </summary>
        public static FluidSimulation CreateOil(AABB bounds)
        {
            return FluidSimulation.Oil(bounds);
        }

        #endregion

        #region Triggers and Sensors

        /// <summary>
        /// Creates a trigger zone.
        /// </summary>
        public static TriggerBody CreateTrigger(Vector3D position, Vector3D halfExtents)
        {
            return TriggerBody.CreateBox(position, halfExtents);
        }

        /// <summary>
        /// Creates a spherical trigger.
        /// </summary>
        public static TriggerBody CreateSphereTrigger(Vector3D position, double radius)
        {
            return TriggerBody.CreateSphere(position, radius);
        }

        /// <summary>
        /// Creates a wind zone sensor.
        /// </summary>
        public static SensorBody CreateWindZone(Vector3D position, Vector3D size, Vector3D windForce)
        {
            return SensorBody.WindZone(position, size, windForce);
        }

        /// <summary>
        /// Creates a slow zone (like water/mud).
        /// </summary>
        public static SensorBody CreateSlowZone(Vector3D position, Vector3D size, double slowFactor = 0.1)
        {
            return SensorBody.SlowZone(position, size, slowFactor);
        }

        #endregion

        #region Smoke Simulation

        /// <summary>
        /// Creates a smoke simulation at a position.
        /// </summary>
        public static SmokeSimulation CreateSmoke(Vector3D position, int maxParticles = 5000)
        {
            return new SmokeSimulation(position, maxParticles);
        }

        /// <summary>
        /// Creates a campfire smoke effect.
        /// </summary>
        public static SmokeSimulation CreateCampfireSmoke(Vector3D position)
        {
            return SmokeSimulation.Campfire(position);
        }

        /// <summary>
        /// Creates a chimney smoke effect.
        /// </summary>
        public static SmokeSimulation CreateChimneySmoke(Vector3D position)
        {
            return SmokeSimulation.Chimney(position);
        }

        /// <summary>
        /// Creates a steam/vapor effect.
        /// </summary>
        public static SmokeSimulation CreateSteam(Vector3D position)
        {
            return SmokeSimulation.Steam(position);
        }

        /// <summary>
        /// Creates an explosion smoke effect.
        /// </summary>
        public static SmokeSimulation CreateExplosionSmoke(Vector3D position, int particleCount = 200)
        {
            return SmokeSimulation.Explosion(position, particleCount);
        }

        #endregion

        #region Fire Simulation

        /// <summary>
        /// Creates a fire simulation at a position.
        /// </summary>
        public static FireSimulation CreateFire(Vector3D position, int maxParticles = 3000)
        {
            return new FireSimulation(position, maxParticles);
        }

        /// <summary>
        /// Creates a campfire.
        /// </summary>
        public static FireSimulation CreateCampfire(Vector3D position)
        {
            return FireSimulation.Campfire(position);
        }

        /// <summary>
        /// Creates a torch fire.
        /// </summary>
        public static FireSimulation CreateTorch(Vector3D position)
        {
            return FireSimulation.Torch(position);
        }

        /// <summary>
        /// Creates a candle flame.
        /// </summary>
        public static FireSimulation CreateCandle(Vector3D position)
        {
            return FireSimulation.Candle(position);
        }

        /// <summary>
        /// Creates a bonfire (large fire).
        /// </summary>
        public static FireSimulation CreateBonfire(Vector3D position)
        {
            return FireSimulation.Bonfire(position);
        }

        /// <summary>
        /// Creates an inferno (intense fire).
        /// </summary>
        public static FireSimulation CreateInferno(Vector3D position)
        {
            return FireSimulation.Inferno(position);
        }

        /// <summary>
        /// Creates a gas burner flame (blue fire).
        /// </summary>
        public static FireSimulation CreateGasBurner(Vector3D position)
        {
            return FireSimulation.GasBurner(position);
        }

        /// <summary>
        /// Creates an explosion fireball.
        /// </summary>
        public static FireSimulation CreateExplosionFire(Vector3D position, double radius = 3.0)
        {
            return FireSimulation.Explosion(position, radius);
        }

        #endregion

        #region Combustion System

        /// <summary>
        /// Creates a combustion system for managing fire spread and extinguishing.
        /// </summary>
        public static CombustionSystem CreateCombustionSystem()
        {
            return new CombustionSystem();
        }

        #endregion

        #region Precomputation and Scientific Computing

        /// <summary>
        /// Creates a precomputed simulation for offline/scientific processing.
        /// </summary>
        public static PrecomputedSimulation CreatePrecomputed(PhysicsWorld world)
        {
            return new PrecomputedSimulation(world);
        }

        /// <summary>
        /// Creates a precomputed simulation with a new world.
        /// </summary>
        public static PrecomputedSimulation CreatePrecomputed()
        {
            return new PrecomputedSimulation();
        }

        /// <summary>
        /// Creates a simulation recorder for data collection.
        /// </summary>
        public static SimulationRecorder CreateRecorder()
        {
            return new SimulationRecorder();
        }

        /// <summary>
        /// Creates a parallel physics processor for high-performance workloads.
        /// </summary>
        /// <param name="threadCount">Number of threads (0 = auto-detect).</param>
        public static ParallelPhysicsProcessor CreateParallelProcessor(int threadCount = 0)
        {
            return new ParallelPhysicsProcessor(threadCount);
        }

        /// <summary>
        /// Creates a Structure of Arrays particle container for SIMD processing.
        /// </summary>
        /// <param name="capacity">Maximum number of particles.</param>
        public static ParticleSoA CreateParticleSoA(int capacity)
        {
            return new ParticleSoA(capacity);
        }

        #endregion

        #region Real-Time Optimized Physics

        /// <summary>
        /// Creates an optimized physics world for real-time applications.
        /// Features: spatial hashing, island solver, warm starting, multi-threading.
        /// </summary>
        /// <param name="spatialCellSize">Cell size for spatial hash (>= largest object diameter).</param>
        public static OptimizedPhysicsWorld CreateOptimizedWorld(double spatialCellSize = 2.0)
        {
            return new OptimizedPhysicsWorld(spatialCellSize);
        }

        /// <summary>
        /// Creates an optimized physics world with custom gravity.
        /// </summary>
        public static OptimizedPhysicsWorld CreateOptimizedWorld(Vector3D gravity, double spatialCellSize = 2.0)
        {
            return new OptimizedPhysicsWorld(gravity, spatialCellSize);
        }

        /// <summary>
        /// Creates a spatial hash for efficient broad-phase collision detection.
        /// </summary>
        public static SpatialHash CreateSpatialHash(double cellSize = 2.0)
        {
            return new SpatialHash(cellSize);
        }

        /// <summary>
        /// Creates a multi-level spatial hash for scenes with varying object sizes.
        /// </summary>
        public static MultiLevelSpatialHash CreateMultiLevelSpatialHash(
            double minSize = 0.5, double maxSize = 16.0, int levels = 4)
        {
            return new MultiLevelSpatialHash(minSize, maxSize, levels);
        }

        /// <summary>
        /// Creates an island solver for parallel constraint solving.
        /// </summary>
        public static IslandSolver CreateIslandSolver(int maxBodies = 10000)
        {
            return new IslandSolver(maxBodies);
        }

        #endregion

        #region GPU Compute

        /// <summary>
        /// Creates a GPU compute accelerator with automatic backend detection.
        /// Supports OpenCL (Intel/AMD/NVIDIA), CUDA (NVIDIA), and CPU fallback.
        /// </summary>
        public static GpuCompute CreateGpuCompute()
        {
            var compute = new GpuCompute();
            compute.Initialize();
            return compute;
        }

        /// <summary>
        /// Creates a GPU compute accelerator with a preferred backend.
        /// </summary>
        public static GpuCompute CreateGpuCompute(GpuBackend preferredBackend)
        {
            var compute = new GpuCompute(preferredBackend);
            compute.Initialize();
            return compute;
        }

        #endregion

        #region Fracture and Destruction

        /// <summary>
        /// Creates a fracture system for realistic object destruction.
        /// </summary>
        public static FractureSystem CreateFractureSystem()
        {
            return new FractureSystem();
        }

        /// <summary>
        /// Creates a fracture system with a specific random seed for reproducible results.
        /// </summary>
        public static FractureSystem CreateFractureSystem(int seed)
        {
            return new FractureSystem(seed);
        }

        /// <summary>
        /// Creates a glass fracture configuration.
        /// </summary>
        public static FractureConfig GlassFracture() => FractureSystem.Glass();

        /// <summary>
        /// Creates a wood fracture configuration.
        /// </summary>
        public static FractureConfig WoodFracture() => FractureSystem.Wood();

        /// <summary>
        /// Creates a stone/concrete fracture configuration.
        /// </summary>
        public static FractureConfig StoneFracture() => FractureSystem.Stone();

        /// <summary>
        /// Creates a brick fracture configuration.
        /// </summary>
        public static FractureConfig BrickFracture() => FractureSystem.Brick();

        /// <summary>
        /// Creates a metal fracture configuration.
        /// </summary>
        public static FractureConfig MetalFracture() => FractureSystem.Metal();

        /// <summary>
        /// Creates an ice fracture configuration.
        /// </summary>
        public static FractureConfig IceFracture() => FractureSystem.Ice();

        /// <summary>
        /// Creates an erodible body (sand, snow, dirt) filling a box shape.
        /// </summary>
        public static ErodibleBody CreateErodibleBox(
            Vector3D center, Vector3D halfExtents, double particleSize = 0.1, double cohesion = 0.5)
        {
            return ErodibleBody.CreateBox(center, halfExtents, particleSize, cohesion);
        }

        /// <summary>
        /// Creates an erodible sphere (sand ball, snowball).
        /// </summary>
        public static ErodibleBody CreateErodibleSphere(
            Vector3D center, double radius, double particleSize = 0.1, double cohesion = 0.5)
        {
            return ErodibleBody.CreateSphere(center, radius, particleSize, cohesion);
        }

        #endregion

        #region Modifiers

        /// <summary>
        /// Creates a modifier system for managing multiple physics modifiers.
        /// </summary>
        public static ModifierSystem CreateModifierSystem()
        {
            return new ModifierSystem();
        }

        /// <summary>
        /// Creates a wind modifier.
        /// </summary>
        public static WindModifier CreateWind(Vector3D direction, double strength = 5.0)
        {
            return new WindModifier(direction, strength);
        }

        /// <summary>
        /// Creates a gentle breeze.
        /// </summary>
        public static WindModifier CreateBreeze(Vector3D direction)
        {
            var wind = WindModifier.Breeze();
            wind.BaseDirection = direction.Normalized();
            return wind;
        }

        /// <summary>
        /// Creates strong wind.
        /// </summary>
        public static WindModifier CreateStrongWind(Vector3D direction)
        {
            var wind = WindModifier.Strong();
            wind.BaseDirection = direction.Normalized();
            return wind;
        }

        /// <summary>
        /// Creates hurricane-force wind.
        /// </summary>
        public static WindModifier CreateHurricane(Vector3D direction)
        {
            var wind = WindModifier.Hurricane();
            wind.BaseDirection = direction.Normalized();
            return wind;
        }

        /// <summary>
        /// Creates a gravity modifier zone.
        /// </summary>
        public static GravityModifier CreateGravityZone(AABB bounds, Vector3D gravity)
        {
            return new GravityModifier { Bounds = bounds, Gravity = gravity };
        }

        /// <summary>
        /// Creates a zero-gravity zone.
        /// </summary>
        public static GravityModifier CreateZeroGravityZone(AABB bounds)
        {
            return GravityModifier.ZeroG(bounds);
        }

        /// <summary>
        /// Creates an attractor point.
        /// </summary>
        public static AttractorModifier CreateAttractor(Vector3D position, double strength = 100.0, double range = 20.0)
        {
            return new AttractorModifier
            {
                Position = position,
                Strength = strength,
                Range = range
            };
        }

        /// <summary>
        /// Creates a repeller point.
        /// </summary>
        public static AttractorModifier CreateRepeller(Vector3D position, double strength = 100.0, double range = 20.0)
        {
            return new AttractorModifier
            {
                Position = position,
                Strength = strength,
                Range = range,
                Repel = true
            };
        }

        /// <summary>
        /// Creates a tornado vortex.
        /// </summary>
        public static VortexModifier CreateTornado(Vector3D position)
        {
            return VortexModifier.Tornado(position);
        }

        /// <summary>
        /// Creates a whirlpool vortex.
        /// </summary>
        public static VortexModifier CreateWhirlpool(Vector3D position)
        {
            return VortexModifier.Whirlpool(position);
        }

        /// <summary>
        /// Creates a turbulence zone.
        /// </summary>
        public static TurbulenceModifier CreateTurbulence(AABB? bounds = null, double strength = 5.0)
        {
            return new TurbulenceModifier { Bounds = bounds, Strength = strength };
        }

        /// <summary>
        /// Creates a water drag zone.
        /// </summary>
        public static DragModifier CreateWaterZone(AABB bounds)
        {
            return DragModifier.Water(bounds);
        }

        /// <summary>
        /// Creates a mud/quicksand zone.
        /// </summary>
        public static DragModifier CreateMudZone(AABB bounds)
        {
            return DragModifier.Mud(bounds);
        }

        #endregion

        #region Object Pools

        /// <summary>
        /// Creates a generic object pool.
        /// </summary>
        public static ObjectPool<T> CreatePool<T>(
            Func<T>? factory = null,
            Action<T>? reset = null,
            int initialSize = 0) where T : class, new()
        {
            return new ObjectPool<T>(factory, reset, initialSize);
        }

        /// <summary>
        /// Creates a list pool.
        /// </summary>
        public static ListPool<T> CreateListPool<T>(int initialCapacity = 16)
        {
            return new ListPool<T>(initialCapacity);
        }

        #endregion
    }
}
