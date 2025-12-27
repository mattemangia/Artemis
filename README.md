# Artemis Physics Engine

**Artemis** is a lightweight, extensible physics engine library for C# applications. Designed for scientific simulations and game development, it provides a clean API for real-time physics similar to NVIDIA PhysX but fully managed in C#.

## Features

- **Cross-Platform**: Compatible with .NET 6/7/8 and .NET Standard 2.1
- **Game Engine Ready**: Works with Unity, MonoGame, Godot (C#), and any C# project
- **Rigid Body Dynamics**: Full 3D physics with collision detection and resolution
- **Real-Time Optimized**: Spatial hashing, island solver, warm starting - like PhysX
- **GPU Compute**: OpenCL (Intel/AMD/NVIDIA), CUDA support with CPU SIMD fallback
- **Particle Systems**: High-performance particle simulation with up to 100,000+ particles
- **Fluid Simulation**: SPH (Smoothed Particle Hydrodynamics) for realistic liquids
- **Smoke & Fire**: Volumetric smoke with turbulence, fire with blackbody radiation colors
- **Combustion System**: Fire spread to flammable objects, water extinguishing
- **Fracture/Destruction**: Realistic shattering with Voronoi, radial, shatter patterns
- **Erodible Objects**: Sand sculptures, snowballs that erode particle-by-particle
- **Cloth Simulation**: Realistic fabric with wind, collisions, and tearing
- **Physics Modifiers**: Wind, gravity zones, attractors, vortex, turbulence
- **Multiple Force Types**: Gravity, wind, magnetism, buoyancy, vortex, explosions, and more
- **Material System**: Elasticity, plasticity, ductility, friction properties
- **SAT Collision**: Separating Axis Theorem for accurate box-box collisions
- **Scientific Computing**: Precomputation, Monte Carlo, parameter sweeps, data export
- **Orbital Mechanics**: Kepler orbits, N-body simulation, Hohmann transfers
- **Propulsion Systems**: Rocket engines, staging, fuel consumption, delta-V
- **G-Force System**: Acceleration tracking, damage accumulation, destruction
- **Planet Gravity**: Game-scale planets (Mario Galaxy style), N-body systems
- **Slope Stability**: Rockfall simulation, weathering, precomputed hazard zones
- **State Transfer**: Bidirectional state capture, time reversal, undo/redo
- **Force Interaction**: Forces combine, amplify, dampen, or cancel each other
- **Ray Tracing**: GPU-accelerated mirror reflections, glass refraction, PBR materials
- **High Performance**: SIMD, multithreading, parallel processing, SoA data layouts
- **Fluent API**: Intuitive, chainable configuration
- **Extensible**: Easy to add custom forces, bodies, and behaviors

---

## Installation

### Option 1: NuGet Package (Recommended)
```bash
dotnet add package Artemis.Physics
```

### Option 2: DLL Reference
1. Download `Artemis.dll` from the [Releases](https://github.com/yourusername/Artemis/releases) page
2. Add the DLL reference to your project:

**Visual Studio:**
- Right-click your project → Add → Reference → Browse
- Select `Artemis.dll`

**Command Line (.csproj):**
```xml
<ItemGroup>
  <Reference Include="Artemis">
    <HintPath>path/to/Artemis.dll</HintPath>
  </Reference>
</ItemGroup>
```

### Option 3: Project Reference
```xml
<ItemGroup>
  <ProjectReference Include="path/to/Artemis.csproj" />
</ItemGroup>
```

---

## Quick Start

### Basic Physics World

```csharp
using Artemis;

// Create a physics world with Earth gravity
var world = Physics.CreateWorld();

// Create a dynamic sphere
var ball = Physics.CreateSphere(
    position: new Vector3D(0, 10, 0),
    radius: 0.5,
    mass: 1.0
);

// Create a static ground
var ground = Physics.CreateStaticBox(
    position: new Vector3D(0, -1, 0),
    halfExtents: new Vector3D(10, 1, 10)
);

// Add bodies to the world
world.AddBody(ball);
world.AddBody(ground);

// Game loop
while (running)
{
    world.Update(deltaTime); // Usually 1/60 for 60 FPS

    // ball.Position now contains updated position
    Console.WriteLine($"Ball position: {ball.Position}");
}
```

---

## Core Concepts

### Bodies

Artemis supports three body types:

| Type | Description | Use Case |
|------|-------------|----------|
| `Dynamic` | Affected by forces and collisions | Players, projectiles, physics objects |
| `Static` | Never moves, infinite mass | Floors, walls, terrain |
| `Kinematic` | User-controlled movement, affects others | Moving platforms, doors |

```csharp
// Dynamic sphere
var ball = Physics.CreateSphere(position, radius: 0.5, mass: 1.0);

// Static box (floor)
var floor = Physics.CreateStaticBox(position, halfExtents);

// Kinematic (controlled manually)
var platform = Physics.CreateBox(position, halfExtents, mass: 0);
platform.BodyType = BodyType.Kinematic;
```

### Materials

Materials define physical properties:

```csharp
// Custom material
var rubber = new PhysicsMaterial("Rubber")
{
    Density = 1100,
    Restitution = 0.85,        // Bounciness (0-1)
    StaticFriction = 1.0,
    DynamicFriction = 0.8,
    YoungsModulus = 1e7,       // Elasticity
    YieldStrength = 1.5e7,     // Plasticity threshold
    Ductility = 0.95           // Deformability (0-1)
};

var ball = Physics.CreateSphere(position, 0.5, 1.0, rubber);

// Or use presets
var metalBall = Physics.CreateSphere(position, 0.5, 5.0, MaterialPresets.Steel());
var woodBlock = Physics.CreateBox(position, halfExtents, 2.0, MaterialPresets.Wood());
```

**Available Material Presets:**
- Metals: `Steel()`, `Aluminum()`, `Copper()`, `Iron()`, `Gold()`
- Non-Metals: `Rubber()`, `Wood()`, `Glass()`, `Concrete()`, `Ice()`
- Granular: `Sand()`, `Gravel()`
- Soft: `SoftTissue()`, `Foam()`
- Game-Friendly: `BouncyBall()`, `Frictionless()`, `Sticky()`

---

## Forces

### Basic Forces

```csharp
// Uniform gravity (default: Earth)
var gravity = Physics.CreateGravity(9.81);
world.AddForce(gravity);

// Custom gravity direction
var sidewaysGravity = new GravityForce(new Vector3D(5, -9.81, 0));
```

### Point Gravity (Planetary)

```csharp
// Attract objects to a point
var planet = Physics.CreatePointGravity(
    position: Vector3D.Zero,
    mass: 1000.0,
    gravitationalConstant: 100.0  // Custom G for game scale
);

// Factory methods
var sun = PointGravityForce.Sun(position);
var earth = PointGravityForce.Earth(position);
var gameGravity = PointGravityForce.GameGravity(position, mass: 100, strength: 2.0);
```

### Wind

```csharp
// Create wind with turbulence
var wind = Physics.CreateWind(
    direction: new Vector3D(1, 0, 0),
    strength: 20.0,
    turbulence: 0.3
);

// Configure gusts using fluent API
wind.WithGusts(frequency: 0.5, strength: 2.0)
    .WithMassScaling(true);

// Presets
var breeze = WindForce.Breeze(Vector3D.Right);
var storm = WindForce.Storm(Vector3D.Right);
```

### Magnetism

```csharp
// Magnet that attracts
var magnet = Physics.CreateMagnet(
    position: new Vector3D(0, 5, 0),
    strength: 100.0,
    attracting: true
);

// Configure as electromagnet
magnet.AsElectromagnet(current: 1.5)
      .WithRange(20.0);

// Presets
var fridge = MagneticForce.Weak(position);
var industrial = MagneticForce.Industrial(position);
```

### Buoyancy (Water/Fluids)

```csharp
// Create buoyancy for water surface at Y=0
var buoyancy = Physics.CreateBuoyancy(
    surfaceHeight: 0,
    fluidDensity: 1000.0  // Water
);

// Add current
buoyancy.WithCurrent(new Vector3D(2, 0, 0))
        .WithDrag(0.5);

// Presets
var saltWater = BuoyancyForce.SaltWater(surfaceHeight: 0);
var mercury = BuoyancyForce.Mercury(surfaceHeight: 0);
```

### Vortex (Whirlpool/Tornado)

```csharp
// Create a vortex
var vortex = Physics.CreateVortex(
    center: new Vector3D(0, 0, 0),
    rotationStrength: 50.0,
    pullStrength: 20.0
);

// Configure
vortex.WithRadius(inner: 0.5, outer: 10.0)
      .WithLift(30.0)
      .WithHeight(20.0)
      .RotatingClockwise(true);

// Presets
var whirlpool = VortexForce.Whirlpool(center, radius: 5);
var tornado = VortexForce.Tornado(center, radius: 10);
var drain = VortexForce.Drain(center);
```

### Explosion

```csharp
// Create explosion
var explosion = Physics.CreateExplosion(
    center: new Vector3D(0, 0, 0),
    force: 1000.0,
    radius: 10.0
);

// Configure
explosion.WithDuration(0.5)
         .WithUpwardBias(0.3)
         .WithShockwave(speed: 50, thickness: 2);

// Trigger
explosion.Explode(new Vector3D(5, 0, 0));

// Update each frame
explosion.Update(deltaTime);

// Presets
var grenade = ExplosionForce.Grenade(center);
var nuclear = ExplosionForce.Nuclear(center);
var implosion = ExplosionForce.Implosion(center);
```

### Spring

```csharp
// Create spring attached to anchor
var spring = Physics.CreateSpring(
    anchor: new Vector3D(0, 10, 0),
    stiffness: 100.0,
    restLength: 2.0,
    damping: 0.1
);

// Bungee (only pulls when stretched)
var bungee = new BungeeForce(anchor, stiffness: 50, restLength: 3);
```

### Custom Force

```csharp
public class MyForce : IForce
{
    public string Id { get; } = "MyForce";
    public bool Enabled { get; set; } = true;

    public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
    {
        // Your force calculation here
        return new Vector3D(0, 100, 0); // Example: constant upward force
    }
}

world.AddForce(new MyForce());
```

---

## Particle Systems

### Basic Particles

```csharp
// Create particle system
var particles = Physics.CreateParticleSystem(maxParticles: 10000);
particles.Gravity = new Vector3D(0, -9.81, 0);
particles.Bounds = new AABB(new Vector3D(-10, 0, -10), new Vector3D(10, 20, 10));

// Spawn particles
particles.Spawn(
    position: new Vector3D(0, 10, 0),
    velocity: new Vector3D(0, -5, 0),
    mass: 0.1,
    radius: 0.05,
    lifetime: 5.0,
    color: 0xFFFF0000  // Red (ARGB)
);

// Spawn burst
particles.SpawnBurst(
    position: new Vector3D(0, 10, 0),
    count: 100,
    spreadAngle: Math.PI,
    minSpeed: 5,
    maxSpeed: 10
);

// Update
particles.Update(deltaTime);

// Access particles for rendering
foreach (var particle in particles.Particles)
{
    if (particle.IsAlive)
    {
        // Render at particle.Position with particle.Color
    }
}
```

### Sand Simulation

```csharp
// Create sand simulation
var bounds = new AABB(
    new Vector3D(-5, 0, -5),
    new Vector3D(5, 10, 5)
);
var sand = Physics.CreateSandSimulation(bounds);

// Add sand ball
sand.AddSandBall(
    center: new Vector3D(0, 8, 0),
    ballRadius: 1.0,
    color: 0xFFFFAA00  // Orange
);

// Update (with substeps for stability)
sand.Update(deltaTime, subSteps: 4);

// Check for connections (for games like Sand Tetris)
bool spansX = sand.CheckSpansAxis(color: 0xFFFFAA00, axis: 0);
```

### Fluid Simulation (SPH)

```csharp
// Create fluid simulation
var fluid = Physics.CreateWater(bounds);

// Add fluid block
fluid.AddFluidBlock(
    min: new Vector3D(-2, 5, -2),
    max: new Vector3D(2, 8, 2),
    spacing: 0.1
);

// Add rigid body colliders
fluid.AddCollider(wall);
fluid.AddCollider(sphere);

// Add external forces
fluid.AddForce(windForce);

// Configure properties
fluid.WithViscosity(0.001)     // Water viscosity
     .WithStiffness(2000)       // How compressible
     .WithSurfaceTension(0.07); // Surface tension

// Update
fluid.Update(deltaTime, subSteps: 2);

// Presets
var water = FluidSimulation.Water(bounds);
var oil = FluidSimulation.Oil(bounds);
var honey = FluidSimulation.Honey(bounds);
var lava = FluidSimulation.Lava(bounds);
```

### Smoke Simulation

```csharp
// Create smoke at position
var smoke = Physics.CreateSmoke(new Vector3D(0, 0, 0));

// Configure
smoke.Buoyancy = 15.0;        // Upward force
smoke.Turbulence = 5.0;       // Random movement
smoke.SpreadRate = 2.0;       // Horizontal expansion
smoke.GrowthRate = 0.5;       // Size increase over time
smoke.Wind = new Vector3D(2, 0, 0);  // Apply wind

// Presets
var campfireSmoke = Physics.CreateCampfireSmoke(position);
var chimneySmoke = Physics.CreateChimneySmoke(position);
var steam = Physics.CreateSteam(position);
var explosionSmoke = Physics.CreateExplosionSmoke(position, particleCount: 200);

// Emit puff manually
smoke.Puff(20);

// Update each frame
smoke.Update(deltaTime);

// Render particles
foreach (var particle in smoke.Particles.Particles)
{
    if (particle.IsAlive)
    {
        // Render smoke particle with fading alpha
    }
}
```

### Fire Simulation

```csharp
// Create fire at position
var fire = Physics.CreateFire(new Vector3D(0, 0, 0));

// Configure
fire.InitialTemperature = 1500;  // Kelvin (affects color)
fire.CoolingRate = 800;          // How fast it cools
fire.Buoyancy = 20;              // Upward force
fire.Turbulence = 8;             // Flickering
fire.FlameHeight = 1.5;          // Height multiplier
fire.Intensity = 1.0;            // 0-1 (can dim/brighten)

// Fire presets
var campfire = Physics.CreateCampfire(position);
var torch = Physics.CreateTorch(position);
var candle = Physics.CreateCandle(position);
var bonfire = Physics.CreateBonfire(position);
var inferno = Physics.CreateInferno(position);
var gasBurner = Physics.CreateGasBurner(position);  // Blue flame

// Explosion fireball
var explosion = Physics.CreateExplosionFire(position, radius: 3.0);

// Fire with linked smoke
var fireWithSmoke = Physics.CreateCampfire(position).WithSmoke();

// Control
fire.Ignite();       // Start burning
fire.Extinguish();   // Stop burning
fire.WithWind(new Vector3D(5, 0, 0));  // Apply wind

// Update each frame
fire.Update(deltaTime);

// Access particles for rendering
foreach (var particle in fire.Particles.Particles)
{
    if (particle.IsAlive)
    {
        // particle.Color contains temperature-based color (blackbody radiation)
        // Hot = white/blue, Medium = yellow/orange, Cool = red
    }
}

// Get fire color at specific temperature
var (r, g, b, a) = FireSimulation.GetFireColor(1200); // Kelvin
```

### Combustion System (Fire Spread & Extinguishing)

```csharp
// Create combustion system
var combustion = Physics.CreateCombustionSystem();

// Register fires
var fire = Physics.CreateCampfire(position);
combustion.RegisterFire(fire);

// Register water (for extinguishing)
var water = Physics.CreateWater(bounds);
combustion.RegisterWater(water);

// Register flammable objects
var woodBox = Physics.CreateBox(position, halfExtents, 5.0, MaterialPresets.Wood());
var combustible = combustion.RegisterCombustible(
    woodBox,
    flammability: 0.7,       // 0-1, how easily catches fire
    burnTime: 30.0,          // Seconds until burned out
    ignitionTemperature: 573 // Kelvin (~300°C)
);

// Or use material presets
combustion.RegisterCombustible(paperStack, FlammableMaterial.Paper);
combustion.RegisterCombustible(oilBarrel, FlammableMaterial.Oil);
combustion.RegisterCombustible(coalPile, FlammableMaterial.Coal);

// Available flammable materials:
// FlammableMaterial.Wood, Paper, Fabric, Oil, Gasoline, Coal, Leaves, Rubber

// Subscribe to events
combustion.OnIgnition += obj => Console.WriteLine($"{obj.Body.Id} caught fire!");
combustion.OnExtinguish += (fire, reason) => Console.WriteLine($"Fire out: {reason}");
combustion.OnBurnedOut += obj => Console.WriteLine($"{obj.Body.Id} burned completely");

// Update each frame
combustion.Update(deltaTime);

// Manual actions
combustion.IgniteObject(combustible);     // Start fire on object
combustion.SplashWater(position, radius: 2.0);  // Extinguish area
combustion.ExtinguishAll();               // Stop all fires

// Queries
var nearest = combustion.GetNearestFire(playerPosition);
bool onFire = combustion.IsPositionOnFire(position, radius: 1.0);
double temp = combustion.GetTemperatureAtPosition(position);  // Kelvin
```

**Fire-Water Interaction:**
- Water particles near fire reduce intensity
- Enough water extinguishes fire completely
- Extinguishing creates steam effect
- Fire cannot spread through water

**Fire Spread:**
- Fire spreads to nearby flammable objects
- Spread rate depends on distance and flammability
- Objects burn down over time
- Burned objects can no longer catch fire

---

## Triggers and Sensors

### Triggers (Detect Overlaps)

```csharp
// Create trigger zone
var trigger = Physics.CreateTrigger(
    position: new Vector3D(0, 1, 0),
    halfExtents: new Vector3D(2, 2, 2)
);

// Subscribe to events
trigger.OnTriggerEnter += body =>
    Console.WriteLine($"{body.Id} entered trigger");

trigger.OnTriggerExit += body =>
    Console.WriteLine($"{body.Id} exited trigger");

trigger.OnTriggerStay += body =>
    Console.WriteLine($"{body.Id} is inside trigger");

// Check overlaps manually
if (trigger.CheckOverlap(player))
{
    // Player is in trigger
}
```

### Sensors (Apply Effects)

```csharp
// Wind zone
var windZone = Physics.CreateWindZone(
    position: new Vector3D(0, 5, 0),
    size: new Vector3D(10, 10, 10),
    windForce: new Vector3D(20, 0, 0)
);

// Slow zone (like water)
var waterZone = Physics.CreateSlowZone(
    position: new Vector3D(0, -2, 0),
    size: new Vector3D(10, 4, 10),
    slowFactor: 0.3
);

// Anti-gravity zone
var floatZone = SensorBody.AntiGravityZone(position, size, strength: 15);

// Apply effects
windZone.ApplyEffect(player);
```

---

## Collision Filtering

```csharp
// Define layers
const int PlayerLayer = 3;
const int EnemyLayer = 4;
const int ProjectileLayer = 5;

// Create body with filter
var player = Physics.CreateSphere(position, 0.5, 1.0);
player.Filter = CollisionFilter.ForLayer(PlayerLayer)
    .WithLayer(EnemyLayer)      // Collides with enemies
    .WithLayer(CollisionLayers.Static);  // Collides with walls

// Projectile that only hits enemies
var bullet = Physics.CreateSphere(position, 0.1, 0.01);
bullet.Filter = CollisionFilter.Create(ProjectileLayer, EnemyLayer);

// Bodies in same group don't collide
player.Filter = player.Filter.InGroup(1);
ally.Filter = ally.Filter.InGroup(1);  // Won't collide with player
```

---

## Collision Events

```csharp
// Subscribe to collision events
world.OnCollision += collision =>
{
    Console.WriteLine($"Collision: {collision.BodyA.Id} hit {collision.BodyB.Id}");
    Console.WriteLine($"  Impact point: {collision.ContactPoint}");
    Console.WriteLine($"  Normal: {collision.Normal}");
    Console.WriteLine($"  Penetration: {collision.Penetration}");

    // Play sound, spawn particles, deal damage, etc.
};
```

---

## Raycasting

```csharp
var hit = world.Raycast(
    origin: cameraPosition,
    direction: cameraForward,
    maxDistance: 100
);

if (hit.HasValue)
{
    Console.WriteLine($"Hit {hit.Value.Body.Id} at {hit.Value.Point}");
    Console.WriteLine($"Distance: {hit.Value.Distance}");
    Console.WriteLine($"Surface normal: {hit.Value.Normal}");
}
```

---

## Unity Integration

```csharp
using UnityEngine;
using Artemis;
using ArtemisVector3D = Artemis.Core.Vector3D;

public class PhysicsManager : MonoBehaviour
{
    private PhysicsWorld _world;
    private Dictionary<string, GameObject> _bodyObjects = new();

    void Start()
    {
        _world = Physics.CreateWorld();

        // Create ground
        var ground = Physics.CreateStaticBox(
            new ArtemisVector3D(0, -0.5, 0),
            new ArtemisVector3D(50, 0.5, 50)
        );
        _world.AddBody(ground);
    }

    public void AddPhysicsObject(GameObject obj, double mass)
    {
        var pos = ToArtemis(obj.transform.position);
        var body = Physics.CreateSphere(pos, 0.5, mass);
        _world.AddBody(body);
        _bodyObjects[body.Id] = obj;
    }

    void FixedUpdate()
    {
        _world.Update(Time.fixedDeltaTime);

        // Sync transforms
        foreach (var body in _world.Bodies)
        {
            if (_bodyObjects.TryGetValue(body.Id, out var obj))
            {
                obj.transform.position = ToUnity(body.Position);
                obj.transform.rotation = ToUnity(body.Rotation);
            }
        }
    }

    // Conversion helpers
    ArtemisVector3D ToArtemis(Vector3 v) => new(v.x, v.y, v.z);
    Vector3 ToUnity(ArtemisVector3D v) => new((float)v.X, (float)v.Y, (float)v.Z);
    UnityEngine.Quaternion ToUnity(Artemis.Core.Quaternion q) =>
        new((float)q.X, (float)q.Y, (float)q.Z, (float)q.W);
}
```

---

## MonoGame Integration

```csharp
using Microsoft.Xna.Framework;
using Artemis;
using ArtemisVector3D = Artemis.Core.Vector3D;

public class Game1 : Game
{
    private PhysicsWorld _physics;

    protected override void Initialize()
    {
        _physics = Physics.CreateWorld();

        // Add physics objects
        var ball = Physics.CreateSphere(new ArtemisVector3D(0, 10, 0), 0.5, 1.0);
        _physics.AddBody(ball);

        base.Initialize();
    }

    protected override void Update(GameTime gameTime)
    {
        double dt = gameTime.ElapsedGameTime.TotalSeconds;
        _physics.Update(dt);

        base.Update(gameTime);
    }

    // Conversion
    Vector3 ToXNA(ArtemisVector3D v) => new((float)v.X, (float)v.Y, (float)v.Z);
}
```

---

## Scientific Computing & Precomputation

Artemis supports offline simulations for scientific analysis, parameter sweeps, and Monte Carlo methods.

### Basic Precomputation

```csharp
// Setup simulation
var world = Physics.CreateWorld();
var ball = Physics.CreateSphere(new Vector3D(0, 100, 0), 0.5, 1.0);
world.AddBody(ball);

// Create precomputed simulation
var sim = Physics.CreatePrecomputed(world);
sim.SubSteps = 4;  // Physics substeps per frame

// Run simulation asynchronously
var result = await sim.RunAsync(
    duration: 10.0,         // Simulate 10 seconds
    timeStep: 1.0 / 60.0    // 60 FPS
);

// Analyze results
Console.WriteLine($"Computed {result.TotalFrames} frames in {result.ComputeTime:F2}s");
Console.WriteLine($"Real-time factor: {result.RealTimeFactor:F1}x");
Console.WriteLine($"Energy conservation error: {result.EnergyConservationError:P2}");

// Access recorded data
var trajectory = result.Recorder.GetTrajectory(ball.Id);
var positions = result.Recorder.GetPositionOverTime(ball.Id);
var velocities = result.Recorder.GetVelocityOverTime(ball.Id);
```

### Parameter Sweep (Parallel)

```csharp
var sim = Physics.CreatePrecomputed();

// Run 100 variations in parallel
var results = await sim.RunParameterSweepAsync(
    worldFactory: i =>
    {
        var world = Physics.CreateWorld();
        var ball = Physics.CreateSphere(
            new Vector3D(0, 10 + i, 0),  // Varying height
            radius: 0.5,
            mass: 1.0 + i * 0.1          // Varying mass
        );
        world.AddBody(ball);
        return world;
    },
    variations: 100,
    duration: 5.0
);

// Analyze all results
foreach (var result in results)
{
    Console.WriteLine($"Variation {result.VariationIndex}: Energy error = {result.EnergyConservationError}");
}
```

### Monte Carlo Simulation

```csharp
var sim = Physics.CreatePrecomputed();

// Run 1000 random samples
var mcResult = await sim.RunMonteCarloAsync(
    worldFactory: random =>
    {
        var world = Physics.CreateWorld();
        var ball = Physics.CreateSphere(
            new Vector3D(
                random.NextDouble() * 10 - 5,  // Random X
                50 + random.NextDouble() * 10, // Random height
                random.NextDouble() * 10 - 5   // Random Z
            ),
            radius: 0.5,
            mass: 1.0
        );
        world.AddBody(ball);
        return world;
    },
    samples: 1000,
    duration: 10.0,
    measurementFunc: result =>
    {
        // Measure final Y position
        var lastFrame = result.Recorder.Frames.Last();
        return lastFrame.BodyStates[0].Position.Y;
    },
    seed: 42  // Reproducible
);

Console.WriteLine($"Mean final Y: {mcResult.Mean:F3} ± {mcResult.StandardDeviation:F3}");
Console.WriteLine($"95% CI: [{mcResult.Mean - mcResult.ConfidenceInterval95:F3}, {mcResult.Mean + mcResult.ConfidenceInterval95:F3}]");
Console.WriteLine($"Range: [{mcResult.Min:F3}, {mcResult.Max:F3}]");
```

### Convergence Analysis

```csharp
var sim = Physics.CreatePrecomputed();

// Test different time steps
var convergence = await sim.RunConvergenceAnalysisAsync(
    worldFactory: () =>
    {
        var world = Physics.CreateWorld();
        world.AddBody(Physics.CreateSphere(new Vector3D(0, 10, 0), 0.5, 1.0));
        return world;
    },
    duration: 5.0,
    timeSteps: new[] { 1.0/30, 1.0/60, 1.0/120, 1.0/240 },
    errorMetric: result => result.EnergyConservationError
);

Console.WriteLine($"Convergence rate: {convergence.ConvergenceRate:F2}");
Console.WriteLine($"Is converging: {convergence.IsConverging}");
```

### High-Performance Parallel Processing

```csharp
// Create parallel processor (auto-detect CPU cores)
var processor = Physics.CreateParallelProcessor();

// Use Structure of Arrays for cache efficiency
var particles = Physics.CreateParticleSoA(capacity: 1_000_000);

// Add particles
for (int i = 0; i < 500_000; i++)
{
    particles.Add(
        position: new Vector3D(random.NextDouble() * 100, 50, random.NextDouble() * 100),
        velocity: Vector3D.Zero,
        mass: 1.0,
        radius: 0.1,
        lifetime: 10.0
    );
}

// Process with SIMD acceleration
processor.UseSIMD = true;
processor.UpdateParticlesSIMD(particles, deltaTime: 0.016, gravity: new Vector3D(0, -9.81, 0));

// Parallel collision detection
var pairs = processor.BroadPhaseParallel(particles, cellSize: 1.0);
processor.ResolveCollisionsParallel(particles, pairs);

// Compact dead particles periodically
particles.Compact();
```

### Data Export

```csharp
var recorder = Physics.CreateRecorder();
recorder.TrackTrajectories = true;
recorder.StartRecording();

// ... run simulation ...

// Export to CSV (for Excel, Python, R, etc.)
recorder.ExportToCsv("simulation_data.csv");
recorder.ExportToCsv("trajectories.csv");
recorder.ExportEnergyToCsv("energy.csv");

// Export to binary (efficient for large datasets)
recorder.ExportToBinary("simulation.artemis");

// Load binary data later
var loaded = SimulationRecorder.ImportFromBinary("simulation.artemis");
```

### Streaming for Large Simulations

```csharp
// For very long simulations, stream directly to disk
var sim = Physics.CreatePrecomputed(world);

await sim.RunStreamingAsync(
    duration: 3600.0,     // 1 hour simulation
    timeStep: 1.0 / 60.0,
    outputPath: "long_simulation.stream"
);

// Read back frame by frame (low memory)
foreach (var frame in PrecomputedSimulation.ReadStreamingData("long_simulation.stream"))
{
    Console.WriteLine($"Frame {frame.FrameIndex}, Time {frame.Time:F3}s");
    foreach (var (bodyId, data) in frame.BodyData)
    {
        Console.WriteLine($"  {bodyId}: pos={data.pos}");
    }
}
```

---

## Real-Time Optimized Physics

For games and real-time applications requiring PhysX-like performance, use `OptimizedPhysicsWorld`:

### Features

- **Spatial Hashing**: O(1) average-case broad-phase collision detection
- **Island Solver**: Groups connected bodies for parallel solving
- **Warm Starting**: Caches impulses for faster convergence
- **Multi-threading**: Parallel collision detection and constraint solving
- **Object Pooling**: Reduces GC pressure for smooth frame rates
- **CCD**: Continuous collision detection for fast-moving objects

### Basic Usage

```csharp
// Create optimized world (cell size should be >= largest object)
var world = Physics.CreateOptimizedWorld(spatialCellSize: 2.0);

// Configure performance settings
world.UseMultiThreading = true;
world.UseWarmStarting = true;
world.UseCCD = true;
world.SolverIterations = 8;       // More = accurate, slower
world.PositionIterations = 3;

// Add bodies as usual
world.AddBody(ball);
world.AddBody(ground);

// Update
world.Update(deltaTime);

// Get performance stats
Console.WriteLine(world.GetPerformanceStats());
// Output:
// Bodies: 1000 (Dynamic: 800)
// Islands: 45 (Sleeping: 12)
// Broad Phase: 0.45ms (2340 pairs)
// Narrow Phase: 1.23ms (156 contacts)
// Solver: 0.89ms
// Integration: 0.21ms
// Spatial Hash Cells: 1234
```

### Spatial Hash for Custom Use

```csharp
// Create spatial hash
var spatialHash = Physics.CreateSpatialHash(cellSize: 2.0);

// Insert bodies
foreach (var body in bodies)
{
    spatialHash.Insert(body);
}

// Query area
var nearby = spatialHash.Query(bounds);

// Get collision pairs (replaces O(n²) with O(n))
var pairs = new List<(IPhysicsBody, IPhysicsBody)>();
spatialHash.GetPotentialPairs(pairs);

// Raycast with spatial acceleration
foreach (var body in spatialHash.RaycastAll(origin, direction, maxDistance))
{
    // Process potential hits
}
```

### Multi-Level Spatial Hash

For scenes with objects of vastly different sizes:

```csharp
// Creates 4 levels from 0.5 to 16.0 cell size
var multiHash = Physics.CreateMultiLevelSpatialHash(
    minSize: 0.5,
    maxSize: 16.0,
    levels: 4
);

// Bodies automatically go to appropriate level based on size
multiHash.Rebuild(bodies);
```

---

## GPU Compute

Accelerate particle physics and SPH fluids with GPU compute. Supports Intel, AMD, and NVIDIA GPUs.

### Automatic Backend Selection

```csharp
// Auto-detect best available backend (CUDA > OpenCL > CPU)
var gpu = Physics.CreateGpuCompute();

Console.WriteLine($"Using: {gpu.ActiveBackend}");
Console.WriteLine($"Device: {gpu.DeviceInfo?.Name}");
Console.WriteLine($"Memory: {gpu.DeviceInfo?.TotalMemory / 1024 / 1024} MB");
```

### Force Specific Backend

```csharp
// Force OpenCL (works on Intel/AMD/NVIDIA)
var gpu = Physics.CreateGpuCompute(GpuBackend.OpenCL);

// Force CUDA (NVIDIA only)
var gpu = Physics.CreateGpuCompute(GpuBackend.CUDA);

// CPU with SIMD (fallback, always available)
var gpu = Physics.CreateGpuCompute(GpuBackend.CPU);
```

### Accelerated Particle Updates

```csharp
var particles = Physics.CreateParticleSoA(100000);
var gpu = Physics.CreateGpuCompute();

// Spawn particles
for (int i = 0; i < 100000; i++)
{
    particles.Add(randomPosition, Vector3D.Zero, 1.0, 0.1, 5.0);
}

// GPU-accelerated update (gravity, integration)
gpu.UpdateParticles(particles, deltaTime: 0.016f, gravity: new Vector3D(0, -9.81, 0));

// GPU-accelerated collision detection
gpu.ComputeParticleCollisions(particles, cellSize: 0.5f);
```

### GPU SPH Fluid Simulation

```csharp
// SPH density and pressure computation on GPU
gpu.ComputeSPH(
    particles,
    smoothingRadius: 0.2f,
    restDensity: 1000f,
    stiffness: 2000f
);
```

### Device Selection

```csharp
var gpu = Physics.CreateGpuCompute();

// List all available devices
foreach (var device in gpu.AvailableDevices)
{
    Console.WriteLine($"[{device.DeviceId}] {device.Name} ({device.Backend})");
    Console.WriteLine($"    Memory: {device.TotalMemory / 1024 / 1024} MB");
    Console.WriteLine($"    Compute Units: {device.ComputeUnits}");
    Console.WriteLine($"    Double Precision: {device.SupportsDouble}");
}

// Select specific device
gpu.SelectDevice(1);
```

---

## Fracture & Destruction

Realistic object shattering with multiple fracture patterns.

### Basic Fracture

```csharp
// Create fracture system
var fracture = Physics.CreateFractureSystem();

// Make objects fracturable
var glassWindow = Physics.CreateBox(position, halfExtents, 10.0, MaterialPresets.Glass());
fracture.MakeFracturable(glassWindow, Physics.GlassFracture());

var woodCrate = Physics.CreateBox(position, halfExtents, 20.0, MaterialPresets.Wood());
fracture.MakeFracturable(woodCrate, Physics.WoodFracture());

// Subscribe to fracture events
fracture.OnFracture += result =>
{
    Console.WriteLine($"Fractured into {result.Fragments.Count} pieces!");

    // Add fragments to physics world
    foreach (var fragment in result.Fragments)
    {
        world.AddBody(fragment);
    }

    // Remove original
    world.RemoveBody(result.OriginalBody);

    // Add debris particles
    foreach (var debris in result.Debris)
    {
        particleSystem.Add(debris);
    }
};

// Check collisions for fracture
world.OnCollision += collision =>
{
    double impulse = CalculateImpulse(collision);
    fracture.CheckAndFracture(collision, impulse);
};
```

### Fracture Patterns

```csharp
// Voronoi (realistic, natural-looking)
var voronoi = new FractureConfig { Pattern = FracturePattern.Voronoi };

// Radial (from impact point outward)
var radial = new FractureConfig { Pattern = FracturePattern.Radial };

// Shatter (many small pieces, glass-like)
var shatter = new FractureConfig { Pattern = FracturePattern.Shatter };

// Splinter (elongated pieces, wood-like)
var splinter = new FractureConfig { Pattern = FracturePattern.Splinter };

// Brick (regular rectangular pieces)
var brick = new FractureConfig { Pattern = FracturePattern.Brick };

// Uniform (grid-based)
var uniform = new FractureConfig { Pattern = FracturePattern.Uniform };
```

### Material Presets

```csharp
var glassConfig = Physics.GlassFracture();   // Many small pieces, low threshold
var woodConfig = Physics.WoodFracture();     // Splinters, cascade fracture
var stoneConfig = Physics.StoneFracture();   // Voronoi chunks
var brickConfig = Physics.BrickFracture();   // Regular pieces
var metalConfig = Physics.MetalFracture();   // Few large pieces, high threshold
var iceConfig = Physics.IceFracture();       // Shatter pattern, very fragile
```

### Custom Fracture Configuration

```csharp
var customConfig = new FractureConfig
{
    Pattern = FracturePattern.Voronoi,
    FractureThreshold = 150.0,         // Minimum impulse to break
    MinPieces = 5,                     // Minimum fragments
    MaxPieces = 20,                    // Maximum fragments
    MinFragmentSize = 0.1,             // Smallest piece (fraction of original)
    GenerateDebris = true,             // Create small debris particles
    DebrisCount = 50,                  // Number of debris particles
    InheritVelocity = true,            // Fragments keep original velocity
    FragmentVelocityScale = 1.5,       // Extra velocity from explosion
    AllowCascadeFracture = true,       // Fragments can further break
    MaxCascadeDepth = 2                // How many times fragments can break
};
```

### Manual Fracture

```csharp
// Trigger fracture manually (e.g., shooting)
var result = fracture.Fracture(
    body: target,
    impactPoint: hitPoint,
    impactDirection: bulletDirection,
    force: 500.0
);

// Process result
foreach (var fragment in result.Fragments)
{
    world.AddBody(fragment);
}
```

---

## Erodible Bodies (Sand, Snow, Dirt)

Create objects made of particles that can erode over time.

### Create Erodible Objects

```csharp
// Sand sculpture (box shape)
var sandCastle = Physics.CreateErodibleBox(
    center: new Vector3D(0, 2, 0),
    halfExtents: new Vector3D(2, 2, 2),
    particleSize: 0.1,
    cohesion: 0.5  // 0 = loose sand, 1 = solid
);

// Snowball (sphere shape)
var snowball = Physics.CreateErodibleSphere(
    center: new Vector3D(5, 1, 0),
    radius: 1.0,
    particleSize: 0.08,
    cohesion: 0.6
);
```

### Wind Erosion (Sand Blown by Wind)

```csharp
// Create wind modifier
var wind = Physics.CreateStrongWind(new Vector3D(1, 0, 0));

// Connect erodible bodies to wind
wind.ErodibleBodies.Add(sandCastle);

// Handle eroded particles
wind.OnParticlesEroded += erodedParticles =>
{
    // Add to particle system for rendering
    foreach (var p in erodedParticles)
    {
        particleSystem.Particles.Add(p);
    }
};

// Update wind each frame
wind.Update(deltaTime);
```

### Impact Erosion (Shooting Sand)

```csharp
// Projectile hits sand sculpture
var eroded = sandCastle.ApplyImpact(
    impactPoint: hitPosition,
    impactForce: 100.0,
    impactRadius: 0.5
);

// Eroded particles fly away
foreach (var p in eroded)
{
    particleSystem.Add(p);
}
```

---

## Physics Modifiers

Higher-level systems that combine forces, particles, and bodies.

### Wind Modifier

```csharp
// Create wind with gusts and turbulence
var wind = Physics.CreateWind(new Vector3D(1, 0, 0), strength: 10.0);

// Or use presets
var breeze = Physics.CreateBreeze(Vector3D.Right);
var strong = Physics.CreateStrongWind(Vector3D.Right);
var hurricane = Physics.CreateHurricane(Vector3D.Right);

// Configure
wind.GustStrength = 2.0;      // Gust multiplier
wind.GustFrequency = 0.5;     // Gusts per second
wind.Turbulence = 0.3;        // Random variation 0-1
wind.Bounds = optionalBounds; // Limit to area

// Apply to bodies and particles
wind.Update(deltaTime);
wind.Apply(body, deltaTime);
wind.Apply(particle, deltaTime);
```

### Gravity Zones

```csharp
// Zero gravity zone
var zeroG = Physics.CreateZeroGravityZone(bounds);

// Custom gravity zone
var moonGravity = Physics.CreateGravityZone(
    bounds: moonArea,
    gravity: new Vector3D(0, -1.62, 0)
);

// Apply
moonGravity.Apply(astronaut, deltaTime);
```

### Attractors & Repellers

```csharp
// Black hole attractor
var blackHole = Physics.CreateAttractor(
    position: new Vector3D(0, 0, 0),
    strength: 500.0,
    range: 30.0
);

// Repeller (force field)
var shield = Physics.CreateRepeller(
    position: playerPosition,
    strength: 200.0,
    range: 5.0
);

// Apply to objects
blackHole.Apply(asteroid, deltaTime);
shield.Apply(enemy, deltaTime);
```

### Vortex Effects

```csharp
// Tornado
var tornado = Physics.CreateTornado(new Vector3D(0, 0, 0));

// Whirlpool (water)
var whirlpool = Physics.CreateWhirlpool(new Vector3D(0, -5, 0));

// Custom vortex
var custom = new VortexModifier
{
    Position = center,
    Axis = Vector3D.Up,
    RotationalStrength = 30.0,
    InwardStrength = 10.0,
    LiftStrength = 20.0,
    OuterRadius = 15.0,
    InnerRadius = 2.0
};
```

### Turbulence

```csharp
var turbulence = Physics.CreateTurbulence(
    bounds: stormArea,
    strength: 8.0
);

turbulence.Frequency = 2.0;  // Oscillation speed
turbulence.Scale = 1.0;      // Spatial scale
```

### Drag Zones (Water, Mud)

```csharp
// Water resistance
var waterZone = Physics.CreateWaterZone(poolBounds);

// Mud/quicksand
var mudZone = Physics.CreateMudZone(swampBounds);
```

---

## G-Force System

Track acceleration and destroy objects under extreme G-forces.

### Basic Usage

```csharp
// Create G-force system
var gForce = Physics.CreateGForceSystem();

// Register bodies with G-force limits
gForce.Register(fighter, Physics.AircraftGForceLimits());
gForce.Register(pilot, Physics.PilotGForceLimits());
gForce.Register(fragilePayload, Physics.FragileGForceLimits());

// Subscribe to destruction events
gForce.OnDestruction += (body, gForces) =>
{
    Console.WriteLine($"{body.Id} destroyed at {gForces}G!");
    // Trigger fracture, explosion, etc.
};

// Update each frame
gForce.Update(deltaTime);
```

### G-Force Limits

```csharp
// Preset limits
var human = Physics.HumanGForceLimits();      // 5G sustained, 9G peak
var pilot = Physics.PilotGForceLimits();      // 9G sustained, 12G peak
var aircraft = Physics.AircraftGForceLimits(); // 9G sustained, 15G peak
var spacecraft = Physics.SpacecraftGForceLimits(); // 3G sustained, 8G peak
var fragile = Physics.FragileGForceLimits();  // 2G sustained, 3G peak

// Custom limits
var custom = new GForceLimits
{
    MaxSustainedG = 8.0f,
    MaxPeakG = 12.0f,
    DamageThresholdG = 6.0f,
    SustainedDuration = 3.0f  // Seconds before damage
};
```

---

## Orbital Mechanics

Simulate celestial bodies, Kepler orbits, and space missions.

### Kepler Orbits

```csharp
// Create orbital mechanics system
var orbital = Physics.CreateOrbitalMechanics();

// Add celestial bodies
var earth = Physics.CreateEarth();
var moon = Physics.CreateMoon();
orbital.AddBody(earth);
orbital.AddBody(moon);

// Calculate orbit parameters
var elements = orbital.CalculateOrbitalElements(satellite, earth);
Console.WriteLine($"Semi-major axis: {elements.SemiMajorAxis}m");
Console.WriteLine($"Eccentricity: {elements.Eccentricity}");
Console.WriteLine($"Period: {elements.Period}s");

// Propagate orbit
var futurePos = orbital.PropagateOrbit(elements, deltaTime);
```

### N-Body Simulation

```csharp
// Create solar system
var sun = Physics.CreateSun();
var earth = Physics.CreateEarth();
var mars = Physics.CreateMars();

orbital.AddBody(sun);
orbital.AddBody(earth);
orbital.AddBody(mars);

// Update (computes gravitational interactions)
orbital.Update(deltaTime);
```

### Hohmann Transfer

```csharp
// Calculate transfer orbit
var transfer = orbital.CalculateHohmannTransfer(
    fromBody: earth,
    toBody: mars,
    spacecraft: probe
);

Console.WriteLine($"Delta-V required: {transfer.TotalDeltaV}m/s");
Console.WriteLine($"Transfer time: {transfer.TransferTime}s");
Console.WriteLine($"Departure burn: {transfer.DepartureDeltaV}m/s");
Console.WriteLine($"Arrival burn: {transfer.ArrivalDeltaV}m/s");
```

---

## Propulsion Systems

Realistic rocket engine simulation with staging and fuel consumption.

### Basic Spacecraft

```csharp
// Create propulsion system
var propulsion = Physics.CreatePropulsionSystem();

// Use preset spacecraft
var falcon9 = Physics.CreateFalcon9();
var saturnV = Physics.CreateSaturnV();
var ionProbe = Physics.CreateIonProbe();

propulsion.AddSpacecraft(falcon9);

// Control throttle
falcon9.CurrentStage.SetThrottle(1.0f);  // Full thrust

// Update
propulsion.Update(deltaTime);

// Check stats
Console.WriteLine($"Delta-V remaining: {falcon9.TotalDeltaV}m/s");
Console.WriteLine($"Fuel: {falcon9.CurrentStage?.PropellantMass}kg");
```

### Custom Spacecraft

```csharp
var spacecraft = new Spacecraft { Name = "Custom Rocket" };

// First stage
var stage1 = new PropulsionStage { Name = "Stage 1", DryMass = 5000 };
stage1.Engines.Add(Engine.Merlin1D());
stage1.Engines.Add(Engine.Merlin1D());
stage1.Propellants.Add(Propellant.RP1(50000));
stage1.Propellants.Add(Propellant.LiquidOxygen(100000));

// Second stage
var stage2 = new PropulsionStage { Name = "Stage 2", DryMass = 1000 };
stage2.Engines.Add(Engine.IonThruster());
stage2.Propellants.Add(Propellant.Xenon(200));

spacecraft.Stages.Add(stage1);
spacecraft.Stages.Add(stage2);
spacecraft.CurrentStageIndex = 0;
stage1.ActivateEngines();
```

### Engine Types

```csharp
var merlin = Physics.CreateMerlin1DEngine();   // Falcon 9 engine
var raptor = Physics.CreateRaptorEngine();      // Starship engine
var ion = Physics.CreateIonThruster();          // High-efficiency ion drive
var rcs = Physics.CreateRCSThruster();          // Attitude control

// Engine presets
Engine.Merlin1D();       // 845 kN, Isp 311s
Engine.Raptor();         // 1850 kN, Isp 380s
Engine.RL10();           // 110 kN, Isp 465s (vacuum)
Engine.IonThruster();    // 0.236 N, Isp 4190s
Engine.HallThruster();   // 0.29 N, Isp 1770s
Engine.SolidBooster();   // 11.8 MN, Isp 268s
```

---

## Planet Gravity

Game-scale planetary gravity (like Mario Galaxy).

### Single Planet

```csharp
// Create game-scale planet
var planet = Physics.CreatePointGravity(
    center: Vector3.Zero,
    strength: 9.81f,
    radius: 10f
);

// Apply to player
Vector3 gravity = planet.GetGravity(playerPosition);
player.Velocity += gravity * deltaTime;

// Get "up" direction (for orientation)
Vector3 up = planet.GetUpDirection(playerPosition);
```

### Multiple Planets

```csharp
// Create multi-planet system
var planets = Physics.CreateMultiPointGravity();

planets.AddPoint(new PointGravity { Center = pos1, Strength = 9.81f, Radius = 5f });
planets.AddPoint(new PointGravity { Center = pos2, Strength = 6.0f, Radius = 3f });
planets.AddPoint(new PointGravity { Center = pos3, Strength = 12.0f, Radius = 8f });

// Gravity uses closest planet
Vector3 gravity = planets.GetGravity(playerPosition, useClosest: true);

// Find which planet player is on
var dominant = planets.GetDominantPoint(playerPosition);
```

### Realistic N-Body

```csharp
// Create N-body system
var gravitySystem = Physics.CreatePlanetGravitySystem();

// Add planets with realistic parameters
gravitySystem.AddBody(Physics.CreateEarthGravity());
gravitySystem.AddBody(Physics.CreateMoonGravity());

// Apply to rigid bodies
gravitySystem.ApplyToRigidBodies(world.Bodies, deltaTime);

// Detect collisions/mergers
var collisions = gravitySystem.DetectCollisions();
foreach (var (a, b) in collisions)
{
    var merged = gravitySystem.MergeBodies(a, b);
}
```

---

## Slope Stability & Rockfall

Simulate landslides, rockfalls, and terrain stability.

### Real-Time Rockfall

```csharp
// Create slope stability system
var slope = Physics.CreateSlopeStabilitySystem();

// Set terrain height function
slope.SetTerrainFunction(
    heightFunc: pos => GetTerrainHeight(pos),
    normalFunc: pos => GetTerrainNormal(pos)
);

// Add slope cells for analysis
for (int x = 0; x < 100; x++)
{
    for (int z = 0; z < 100; z++)
    {
        slope.AddCell(new SlopeCell
        {
            Position = new Vector3(x, GetHeight(x, z), z),
            Normal = GetNormal(x, z),
            Height = 5.0f,
            Mass = 1000f,
            Material = Physics.RockMaterial()
        });
    }
}

// Analyze stability
slope.AnalyzeAllCells();
var unstable = slope.FindUnstableCells(threshold: 1.0f);

// Trigger collapse
foreach (var cell in unstable)
{
    slope.TriggerCollapse(cell, fragmentCount: 15);
}

// Update rockfall
slope.Update(deltaTime);

// Render fragments
foreach (var fragment in slope.ActiveFragments)
{
    if (fragment.IsActive)
        RenderRock(fragment.Position, fragment.Radius);
}
```

### Precomputed Rockfall

```csharp
// Precompute rockfall for hazard mapping
var precomputed = slope.PrecomputeRockfall(
    sourcePosition: cliffTop,
    sourceMass: 5000f,
    material: Physics.RockMaterial(),
    duration: 30f,
    timeStep: 0.02f,
    fragmentCount: 50
);

Console.WriteLine($"Max runout: {precomputed.FinalResult.MaxRunoutDistance}m");
Console.WriteLine($"Max energy: {precomputed.FinalResult.MaxEnergy}J");

// Play forward
for (float t = 0; t < precomputed.TotalDuration; t += 0.1f)
{
    var frame = precomputed.GetFrame(t);
    // Render frame
}

// Play in reverse
for (float t = precomputed.TotalDuration; t >= 0; t -= 0.1f)
{
    var frame = precomputed.GetFrameReverse(t);
    // Render reversed
}
```

### Slope Materials

```csharp
var rock = Physics.RockMaterial();     // Angle of repose: 45°
var gravel = Physics.GravelMaterial(); // Angle of repose: 35°
var sand = Physics.SandMaterial();     // Angle of repose: 34°
var clay = Physics.ClayMaterial();     // Angle of repose: 25°
var soil = Physics.SoilMaterial();     // Angle of repose: 30°
var ice = Physics.IceMaterial();       // Angle of repose: 35°
var snow = Physics.SnowMaterial();     // Angle of repose: 38°
```

### Weather Effects

```csharp
// Apply rainfall (increases saturation, reduces stability)
slope.ApplyRainfall(intensity: 0.1f, deltaTime);

// Apply drainage
slope.ApplyDrainage(rate: 0.02f, deltaTime);

// Apply weathering over time
slope.ApplyWeathering(deltaTime);
```

---

## Cloth Simulation

Realistic cloth simulation with force interaction, collisions, and tearing.

### Basic Usage

```csharp
// Create cloth mesh (40x30 particles, 0.1m spacing)
var cloth = Physics.CreateCloth(
    origin: new Vector3(0, 5, 0),
    width: 40,
    height: 30,
    spacing: 0.1f,
    material: Physics.CottonCloth()
);

// Pin top edge (like hanging from a rod)
cloth.PinTopEdge();

// Or pin just corners
cloth.PinCorners();

// Set gravity and wind
cloth.Gravity = new Vector3(0, -9.81f, 0);
cloth.Wind = new Vector3(5, 0, 2);
cloth.WindTurbulence = 0.3f;

// Update each frame
cloth.Update(deltaTime);

// Render vertices
foreach (var (pos, normal, uv) in cloth.GetRenderVertices())
{
    // Draw triangle with pos, normal, uv
}
```

### Cloth Materials

```csharp
var cotton = Physics.CottonCloth();     // Light, flexible
var silk = Physics.SilkCloth();         // Very light, smooth
var denim = Physics.DenimCloth();       // Heavy, stiff
var leather = Physics.LeatherCloth();   // Very stiff, durable
var rubber = Physics.RubberCloth();     // Stretchy
var flag = Physics.FlagCloth();         // Light, no tearing
var sail = Physics.SailCloth();         // Strong, high wind resistance
var curtain = Physics.CurtainCloth();   // Medium weight
```

### Force Interaction

```csharp
// Apply wind from a wind modifier
var wind = Physics.CreateStrongWind(Vector3D.Right);
cloth.ApplyWind(wind);

// Add custom force callback
cloth.AddForce(position => {
    // Example: radial force from explosion
    Vector3 toExplosion = explosionCenter - position;
    float dist = toExplosion.Length();
    if (dist < explosionRadius)
    {
        return Vector3.Normalize(toExplosion) * explosionForce * (1 - dist / explosionRadius);
    }
    return Vector3.Zero;
});

// Apply impulse at a point (e.g., from projectile)
cloth.ApplyImpulse(hitPoint, impulseVector, radius: 0.5f);
```

### Collisions

```csharp
// Add sphere collider (e.g., character's head)
cloth.AddSphereCollider(center: headPosition, radius: 0.15f);

// Moving sphere (update position each frame)
cloth.AddSphereCollider(Vector3.Zero, 0.2f, _ => characterPosition);

// Add plane collider (e.g., floor)
cloth.AddPlaneCollider(point: Vector3.Zero, normal: Vector3.UnitY);

// Add box collider
cloth.AddBoxCollider(min: new Vector3(-1, 0, -1), max: new Vector3(1, 2, 1));

// Handle collision events
cloth.OnCollision += collision =>
{
    Console.WriteLine($"Cloth hit at {collision.ContactPoint}");
};
```

### Tearing

```csharp
// Enable tearing
cloth.Material.CanTear = true;
cloth.Material.TearThreshold = 5f;  // Lower = easier to tear

// Manual tear along a line
cloth.Tear(
    start: new Vector3(0, 3, 0),
    end: new Vector3(2, 3, 0),
    width: 0.1f
);

// Handle tear events
cloth.OnConstraintBreak += constraint =>
{
    Console.WriteLine("Cloth torn!");
    // Spawn particles, play sound, etc.
};
```

### Interactive Cloth

```csharp
// Move pinned particle (e.g., for dragging)
cloth.MovePinnedParticle(x: 0, y: 0, newPosition: mouseWorldPosition);

// Pin/unpin dynamically
cloth.PinParticle(x, y);
cloth.UnpinParticle(x, y);
```

---

## State Transfer & Time Reversal

Capture, save, and reverse physics states.

### Basic State Capture

```csharp
// Create state transfer system
var stateTransfer = Physics.CreateStateTransfer();

// Capture current state
var state = stateTransfer.CaptureState(world.Bodies, currentTime, "Checkpoint1");

// Later, restore state
stateTransfer.ApplyState(state, world.Bodies);
```

### Undo/Redo

```csharp
// Push states to history
stateTransfer.PushState(stateTransfer.CaptureState(bodies, time));

// Undo (go back)
if (stateTransfer.CanUndo)
{
    var prevState = stateTransfer.Undo();
    stateTransfer.ApplyState(prevState, bodies);
}

// Redo (go forward)
if (stateTransfer.CanRedo)
{
    var nextState = stateTransfer.Redo();
    stateTransfer.ApplyState(nextState, bodies);
}
```

### Time Reversal

```csharp
// Capture final state
var finalState = stateTransfer.CaptureState(bodies, time);
finalState.Metadata.IsFinalState = true;

// Create reversed state (velocities negated)
var reversed = stateTransfer.CreateReversedState(finalState);

// Use final as initial (for playing simulation backwards)
var initialFromFinal = stateTransfer.FinalToInitial(finalState, reverseVelocities: true);

// Apply and run simulation backwards
stateTransfer.ApplyState(initialFromFinal, bodies);
```

### Recording and Playback

```csharp
// Create recorder
var recorder = Physics.CreateSimulationRecorder();
recorder.RecordInterval = 0.02f;  // 50 FPS

// Start recording
recorder.StartRecording();

// During simulation
recorder.RecordFrame(stateTransfer.CaptureState(bodies, time), deltaTime);

// Stop recording
recorder.StopRecording();

// Playback forward
foreach (var state in recorder.PlayForward(recorder.GetAllFrames()))
{
    stateTransfer.ApplyState(state, bodies);
    Render();
}

// Playback reverse
foreach (var state in recorder.PlaybackReverse(stateTransfer))
{
    stateTransfer.ApplyState(state, bodies);
    Render();
}
```

### Save/Load State

```csharp
// Save state to file
stateTransfer.SaveToFile(state, "checkpoint.state", compress: true);

// Load state from file
var loaded = stateTransfer.LoadFromFile("checkpoint.state");
stateTransfer.ApplyState(loaded, bodies);
```

---

## Force Interaction

Forces can interact with each other - combining, amplifying, dampening, or canceling.

### Composite Forces

```csharp
// Create composite force that combines multiple forces
var composite = Physics.CreateCompositeForce();

// Add forces
composite.Add(Physics.CreateGravity());
composite.Add(Physics.CreateWind(Vector3D.Right, 10));
composite.Add(Physics.CreateMagnet(center, 100));

// Set interaction rules
composite.SetRule<WindForce, BuoyancyForce>(ForceBlendMode.Additive);
composite.SetAmplify<MagneticForce>(factor: 1.5);  // Same-type forces amplify
composite.SetDampen<GravityForce, BuoyancyForce>(factor: 0.5);  // Opposing forces dampen

// Calculate combined force
Vector3D totalForce = composite.Calculate(position, velocity, mass);
```

### Force Interaction System

```csharp
// Create force interaction system
var forces = Physics.CreateForceInteractionSystem();

// Add global forces
forces.AddGlobalForce(Physics.CreateGravity());

// Create force groups
var windGroup = forces.CreateGroup("wind");
windGroup.Add(Physics.CreateWind(Vector3D.Right, 5));
windGroup.Add(Physics.CreateWind(Vector3D.Forward, 3));

var magneticGroup = forces.CreateGroup("magnetic");
magneticGroup.Add(Physics.CreateMagnet(pos1, 100));
magneticGroup.Add(Physics.CreateMagnet(pos2, 50));

// Enable interactions
forces.EnableWindFireInteraction(amplification: 2.0);  // Wind fans flames
forces.EnableMagneticInterference();  // Opposing magnets cancel
forces.EnableMultiGravity();  // Multiple gravity sources combine

// Calculate total force at position
Vector3D total = forces.CalculateTotalForce(position, velocity, mass);
```

### Blend Modes

```csharp
// Available blend modes
ForceBlendMode.Additive     // Add forces together (default)
ForceBlendMode.Maximum      // Use strongest force only
ForceBlendMode.Minimum      // Use weakest force only
ForceBlendMode.Average      // Average all forces
ForceBlendMode.Multiply     // Multiply forces (for scaling)
ForceBlendMode.Interference // Opposing forces cancel each other
ForceBlendMode.Priority     // Higher priority forces override
```

### Presets

```csharp
// Weather system (wind + thermal interaction)
var weather = Physics.CreateWeatherForces();
weather.GetOrCreateGroup("wind").Add(Physics.CreateStrongWind(Vector3D.Right));

// Space (multi-gravity bodies)
var space = Physics.CreateSpaceForces();
space.AddGlobalForce(Physics.CreatePointGravity(sunPos, sunMass));
space.AddGlobalForce(Physics.CreatePointGravity(planetPos, planetMass));

// Fluid (drag + buoyancy)
var fluid = Physics.CreateFluidForces();

// Electromagnetic
var em = Physics.CreateElectromagneticForces();
```

---

## Ray Tracing

GPU-accelerated ray tracing for mirror reflections, glass refraction, and realistic materials.

### Basic Usage

```csharp
// Create ray tracing system (auto-detects GPU)
var rt = Physics.CreateRayTracing();

// Add materials
int mirrorMat = rt.AddMaterial(Physics.MirrorMaterial());
int glassMat = rt.AddMaterial(Physics.GlassMaterial());
int goldMat = rt.AddMaterial(Physics.GoldMaterial());

// Add objects
rt.AddObject(new TraceableSphere
{
    Center = new Vector3(0, 1, 0),
    Radius = 1.0f,
    MaterialIndex = mirrorMat
});

rt.AddObject(new TraceablePlane
{
    Center = Vector3.Zero,
    Normal = Vector3.UnitY,
    HalfExtents = new Vector2(10, 10),
    MaterialIndex = goldMat
});

// Trace single ray
var color = rt.TraceRay(new Ray(origin, direction));

// Trace full image
var pixels = rt.TraceImage(
    width: 1920,
    height: 1080,
    cameraPos: new Vector3(0, 2, -5),
    cameraTarget: Vector3.Zero,
    cameraUp: Vector3.UnitY,
    fov: 60f
);

// Convert to byte array for display
byte[] imageData = rt.ToByteArray(pixels, exposure: 1.0f);
```

### Mirror Objects

```csharp
// Add mirror sphere
rt.AddMirrorSphere(center: new Vector3(0, 1, 0), radius: 1.0f);

// Add mirror plane (wall mirror)
rt.AddMirrorPlane(
    center: new Vector3(0, 2, 5),
    normal: new Vector3(0, 0, -1),
    size: new Vector2(4, 3)
);

// Add mirror box (reflective cube)
rt.AddMirrorBox(center: new Vector3(2, 0.5f, 0), size: new Vector3(1, 1, 1));
```

### Materials

```csharp
// Preset materials
var mirror = Physics.MirrorMaterial();    // Perfect reflection
var chrome = Physics.ChromeMaterial();    // Polished metal
var gold = Physics.GoldMaterial();        // Gold metal
var copper = Physics.CopperMaterial();    // Copper metal
var glass = Physics.GlassMaterial();      // Transparent glass
var water = Physics.WaterMaterial();      // Water surface

// Custom diffuse material
var redMatte = Physics.DiffuseMaterial(new Vector3(1, 0, 0));

// Custom glossy material
var blueGlossy = Physics.GlossyMaterial(new Vector3(0, 0, 1), roughness: 0.3f);

// Fully custom material
var custom = new RayTracingMaterial
{
    Albedo = new Vector3(0.8f, 0.2f, 0.1f),  // Color
    Reflectivity = 0.7f,                       // 0 = diffuse, 1 = mirror
    Roughness = 0.1f,                          // 0 = smooth, 1 = rough
    Metalness = 0.5f,                          // 0 = dielectric, 1 = metal
    Transparency = 0.0f,                       // 0 = opaque, 1 = transparent
    IOR = 1.5f                                 // Index of refraction
};
```

### Lights

```csharp
// Point light
var point = Physics.CreatePointLight(
    position: new Vector3(5, 10, 5),
    color: new Vector3(1, 1, 1)
);

// Directional light (sun)
var sun = Physics.CreateDirectionalLight(
    direction: new Vector3(-1, -1, -1),
    color: new Vector3(1, 0.95f, 0.8f)
);

// Area light (soft shadows)
var area = Physics.CreateAreaLight(
    position: new Vector3(0, 5, 0),
    radius: 2.0f,
    color: new Vector3(1, 1, 1)
);
```

### GPU Backends

```csharp
// Auto-detect best backend
var rt = Physics.CreateRayTracing();

// Force specific backend
var rtOpenCL = Physics.CreateRayTracing(GpuBackend.OpenCL);  // Intel/AMD/NVIDIA
var rtCUDA = Physics.CreateRayTracing(GpuBackend.CUDA);      // NVIDIA only

// CPU-only (always works)
var rtCPU = Physics.CreateRayTracingCPU();

// Check device
Console.WriteLine($"Using: {rt.GpuDevice?.Name}");
```

### Configuration

```csharp
var rt = Physics.CreateRayTracing();

// Settings
rt.MaxBounces = 8;                                    // Reflection depth
rt.BackgroundColor = new Vector3(0.1f, 0.1f, 0.15f); // Sky color
rt.AmbientLight = new Vector3(0.1f, 0.1f, 0.1f);     // Ambient lighting
rt.UseGPU = true;                                     // GPU acceleration
```

---

### Modifier System

Manage multiple modifiers together:

```csharp
var modifiers = Physics.CreateModifierSystem();

// Add modifiers
modifiers.Add(wind);
modifiers.Add(tornado);
modifiers.Add(zeroGravity);

// Update all modifiers
modifiers.Update(deltaTime);

// Apply to all bodies
modifiers.ApplyTo(world.Bodies, deltaTime);

// Apply to particle system
modifiers.ApplyTo(particleSystem, deltaTime);
```

---

## Example: Sand Sculpture Blown by Wind

```csharp
// Create sand sculpture
var sandSculpture = Physics.CreateErodibleBox(
    center: new Vector3D(0, 1, 0),
    halfExtents: new Vector3D(1, 1, 1),
    particleSize: 0.05,
    cohesion: 0.4  // Somewhat loose
);

// Create particle system for falling sand
var fallenSand = Physics.CreateParticleSystem(50000);
fallenSand.Gravity = new Vector3D(0, -9.81, 0);

// Create wind
var wind = Physics.CreateStrongWind(new Vector3D(1, 0, 0));
wind.ErodibleBodies.Add(sandSculpture);

// Handle eroded particles
wind.OnParticlesEroded += particles =>
{
    foreach (var p in particles)
    {
        fallenSand.Spawn(p.Position, p.Velocity, p.Mass, 0.05, p.Lifetime, 0xFFE0C080);
    }
};

// Game loop
while (running)
{
    wind.Update(deltaTime);
    fallenSand.Update(deltaTime);

    // Render sand sculpture particles
    foreach (var p in sandSculpture.Particles)
    {
        DrawParticle(p.Position, p.Color);
    }

    // Render fallen sand
    foreach (var p in fallenSand.Particles)
    {
        if (p.IsAlive)
            DrawParticle(p.Position, 0xFFE0C080);
    }
}
```

---

## Example: Shooting Objects to Break Them

```csharp
// Setup fracture system
var fracture = Physics.CreateFractureSystem();

// Create targets
var targets = new List<RigidBody>();
for (int i = 0; i < 10; i++)
{
    var target = Physics.CreateBox(
        new Vector3D(i * 3, 2, 10),
        new Vector3D(0.5, 0.5, 0.5),
        5.0,
        MaterialPresets.Wood()
    );
    world.AddBody(target);
    fracture.MakeFracturable(target, Physics.WoodFracture());
    targets.Add(target);
}

// Handle fracture
fracture.OnFracture += result =>
{
    world.RemoveBody(result.OriginalBody);
    targets.Remove(result.OriginalBody as RigidBody);

    foreach (var fragment in result.Fragments)
    {
        world.AddBody(fragment);
    }

    // Spawn debris particles
    foreach (var debris in result.Debris)
    {
        particleSystem.Add(debris);
    }
};

// Shooting (e.g., on mouse click)
void Shoot(Vector3D origin, Vector3D direction)
{
    var hit = world.Raycast(origin, direction, 100);
    if (hit.HasValue && hit.Value.Body is RigidBody rb)
    {
        // Apply impact force
        double impactForce = 300.0;
        rb.ApplyImpulseAtPoint(direction * impactForce, hit.Value.Point);

        // Check for fracture
        fracture.Fracture(rb, hit.Value.Point, direction, impactForce);
    }
}
```

---

## Performance Tips

1. **Use spatial partitioning**: For >100 bodies, the engine uses spatial hashing automatically
2. **Enable sleeping**: `body.CanSleep = true` for resting objects
3. **Use layers**: Filter collisions to reduce checks
4. **Substeps**: More substeps = more stable but slower
5. **Particle limits**: Keep particle counts reasonable (<50,000 for real-time)
6. **Fluid resolution**: Larger smoothing radius = fewer particles needed
7. **Use ParticleSoA**: For >10,000 particles, use Structure of Arrays layout
8. **Enable SIMD**: `processor.UseSIMD = true` for vectorized math
9. **Parallel processing**: Use `CreateParallelProcessor()` for multi-core systems
10. **Streaming export**: Use `RunStreamingAsync()` for long simulations to avoid memory issues

---

## API Reference

### Static Factory Methods (Physics class)

| Method | Description |
|--------|-------------|
| **Core** | |
| `CreateWorld()` | Create physics world with Earth gravity |
| `CreateOptimizedWorld(cellSize)` | Create optimized world with spatial hashing |
| `CreateSphere(pos, radius, mass)` | Create dynamic sphere |
| `CreateBox(pos, halfExtents, mass)` | Create dynamic box |
| `CreateStaticBox(pos, halfExtents)` | Create static box |
| **Forces** | |
| `CreateGravity(magnitude)` | Create gravity force |
| `CreateWind(dir, strength, turbulence)` | Create wind force |
| `CreateMagnet(pos, strength, attracting)` | Create magnetic force |
| `CreateBuoyancy(surface, density)` | Create buoyancy force |
| `CreateVortex(center, rotation, pull)` | Create vortex force |
| `CreateExplosion(center, force, radius)` | Create explosion force |
| **Particles & Fluids** | |
| `CreateParticleSystem(max)` | Create particle system |
| `CreateFluidSimulation(bounds)` | Create SPH fluid |
| `CreateSmoke(pos)` | Create smoke simulation |
| `CreateFire(pos)` | Create fire simulation |
| `CreateCampfire(pos)` | Create campfire preset |
| `CreateCombustionSystem()` | Create fire/water interaction system |
| **Real-Time Optimizations** | |
| `CreateSpatialHash(cellSize)` | Create spatial hash for collision detection |
| `CreateMultiLevelSpatialHash()` | Multi-level spatial hash for varied sizes |
| `CreateIslandSolver(maxBodies)` | Create island solver for parallel solving |
| **GPU Compute** | |
| `CreateGpuCompute()` | Create GPU accelerator (auto-detect backend) |
| `CreateGpuCompute(backend)` | Create GPU with specific backend |
| **Destruction** | |
| `CreateFractureSystem()` | Create fracture system for shattering |
| `GlassFracture()` | Glass fracture configuration |
| `WoodFracture()` | Wood fracture configuration |
| `StoneFracture()` | Stone/concrete fracture configuration |
| `MetalFracture()` | Metal fracture configuration |
| `IceFracture()` | Ice fracture configuration |
| `CreateErodibleBox()` | Create erodible sand/snow box |
| `CreateErodibleSphere()` | Create erodible sphere |
| **Modifiers** | |
| `CreateModifierSystem()` | Create modifier manager |
| `CreateWind(dir, strength)` | Create wind modifier |
| `CreateBreeze(dir)` | Create gentle breeze |
| `CreateStrongWind(dir)` | Create strong wind |
| `CreateHurricane(dir)` | Create hurricane-force wind |
| `CreateGravityZone(bounds, gravity)` | Create gravity zone |
| `CreateZeroGravityZone(bounds)` | Create zero-G zone |
| `CreateAttractor(pos, strength, range)` | Create attractor point |
| `CreateRepeller(pos, strength, range)` | Create repeller point |
| `CreateTornado(pos)` | Create tornado vortex |
| `CreateWhirlpool(pos)` | Create whirlpool vortex |
| `CreateTurbulence(bounds, strength)` | Create turbulence zone |
| `CreateWaterZone(bounds)` | Create water drag zone |
| `CreateMudZone(bounds)` | Create mud/quicksand zone |
| **Scientific Computing** | |
| `CreatePrecomputed(world)` | Create precomputed simulation |
| `CreateRecorder()` | Create simulation recorder |
| `CreateParallelProcessor()` | Create parallel physics processor |
| `CreateParticleSoA(capacity)` | Create SoA particle container |
| **G-Force System** | |
| `CreateGForceSystem()` | Create G-force tracking system |
| `HumanGForceLimits()` | G-force limits for humans (5G sustained) |
| `PilotGForceLimits()` | G-force limits for pilots (9G sustained) |
| `AircraftGForceLimits()` | G-force limits for aircraft (9G sustained) |
| `SpacecraftGForceLimits()` | G-force limits for spacecraft |
| `FragileGForceLimits()` | G-force limits for fragile objects |
| **Orbital Mechanics** | |
| `CreateOrbitalMechanics()` | Create orbital mechanics simulation |
| `CreateCelestialBody(name, mass, radius, pos)` | Create celestial body |
| `CreateEarth()` | Create Earth celestial body |
| `CreateMoon()` | Create Moon celestial body |
| `CreateSun()` | Create Sun celestial body |
| `CreateMars()` | Create Mars celestial body |
| **Propulsion Systems** | |
| `CreatePropulsionSystem()` | Create propulsion system manager |
| `CreateFalcon9()` | Create Falcon 9 rocket preset |
| `CreateSaturnV()` | Create Saturn V rocket preset |
| `CreateIonProbe()` | Create ion propulsion probe |
| `CreateMerlin1DEngine()` | Create Merlin 1D engine |
| `CreateRaptorEngine()` | Create Raptor engine |
| `CreateIonThruster()` | Create ion thruster |
| `CreateRCSThruster()` | Create RCS thruster |
| **Planet Gravity** | |
| `CreatePlanetGravitySystem()` | Create N-body gravity system |
| `CreateGravitationalBody(mass, radius, pos)` | Create gravitational body |
| `CreateEarthGravity()` | Create Earth gravity body |
| `CreateMoonGravity()` | Create Moon gravity body |
| `CreateGamePlanet(radius, surfaceGravity)` | Create game-scale planet |
| `CreatePointGravity(center, strength, radius)` | Create point gravity (Mario Galaxy style) |
| `CreateMultiPointGravity()` | Create multi-point gravity system |
| **Slope Stability** | |
| `CreateSlopeStabilitySystem()` | Create slope stability system |
| `RockMaterial()` | Rock slope material (45° repose) |
| `GravelMaterial()` | Gravel slope material (35° repose) |
| `SandMaterial()` | Sand slope material (34° repose) |
| `ClayMaterial()` | Clay slope material (25° repose) |
| `SoilMaterial()` | Soil slope material (30° repose) |
| `IceMaterial()` | Ice slope material (35° repose) |
| `SnowMaterial()` | Snow slope material (38° repose) |
| **State Transfer** | |
| `CreateStateTransfer()` | Create state transfer system |
| `CreateSimulationRecorder()` | Create simulation recorder for playback |
| **Cloth Simulation** | |
| `CreateCloth(origin, width, height, spacing)` | Create cloth mesh |
| `CottonCloth()` | Cotton cloth material |
| `SilkCloth()` | Silk cloth material |
| `DenimCloth()` | Denim cloth material |
| `LeatherCloth()` | Leather material |
| `RubberCloth()` | Rubber material |
| `FlagCloth()` | Flag material (no tearing) |
| `SailCloth()` | Sail material (high wind resistance) |
| `CurtainCloth()` | Curtain material |
| **Force Interaction** | |
| `CreateForceInteractionSystem()` | Create force interaction manager |
| `CreateCompositeForce()` | Create composite force with interaction rules |
| `CreateWeatherForces()` | Weather force preset (wind + thermal) |
| `CreateSpaceForces()` | Space force preset (multi-gravity) |
| `CreateFluidForces()` | Fluid force preset (drag + buoyancy) |
| `CreateElectromagneticForces()` | Electromagnetic force preset |
| **Ray Tracing** | |
| `CreateRayTracing()` | Create GPU ray tracing system |
| `CreateRayTracing(backend)` | Create ray tracing with specific backend |
| `CreateRayTracingCPU()` | Create CPU-only ray tracing |
| `MirrorMaterial()` | Perfect mirror material |
| `ChromeMaterial()` | Chrome/polished metal material |
| `GoldMaterial()` | Gold material |
| `CopperMaterial()` | Copper material |
| `GlassMaterial()` | Glass material (refraction) |
| `WaterMaterial()` | Water material (refraction) |
| `DiffuseMaterial(color)` | Diffuse/matte material |
| `GlossyMaterial(color, roughness)` | Glossy material |
| `CreateTraceableSphere()` | Create traceable sphere |
| `CreateTraceablePlane()` | Create traceable plane |
| `CreateTraceableBox()` | Create traceable box |
| `CreatePointLight()` | Create point light |
| `CreateDirectionalLight()` | Create directional light |
| `CreateAreaLight()` | Create area light (soft shadows) |
| **Object Pools** | |
| `CreatePool<T>()` | Create generic object pool |
| `CreateListPool<T>()` | Create list pool |

---

## License

MIT License - Free for commercial and non-commercial use.

---

## Contributing

Contributions welcome! Please read our contributing guidelines and submit PRs.

---

## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/Artemis/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/Artemis/discussions)
