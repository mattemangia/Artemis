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
    }
}
