using ArtemisEngine;
using AdvancedPhysicsDemo;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== Advanced Physics Demo ===");
        Console.WriteLine("Artemis Physics Engine - Advanced Features Showcase");
        Console.WriteLine();
        Console.WriteLine("This demo showcases:");
        Console.WriteLine("  - Joints & Constraints (chains, vehicles, ragdolls)");
        Console.WriteLine("  - Raycasting");
        Console.WriteLine("  - Trigger Zones");
        Console.WriteLine("  - Collision Events & Layers");
        Console.WriteLine("  - Spatial Partitioning & Sleeping");
        Console.WriteLine();
        Console.WriteLine("Select a demo:");
        Console.WriteLine("  1 - Ragdoll Physics");
        Console.WriteLine("  2 - Vehicle with Suspension");
        Console.WriteLine("  3 - Chain & Bridge");
        Console.WriteLine("  4 - Raycast Shooting");
        Console.WriteLine("  5 - Trigger Zones & Events");
        Console.WriteLine();
        Console.Write("Choice (1-5): ");

        var choice = Console.ReadKey().KeyChar;
        Console.WriteLine();
        Console.WriteLine();

        Console.CursorVisible = false;
        Console.Clear();

        Game? game = choice switch
        {
            '1' => new RagdollDemo(),
            '2' => new VehicleDemo(),
            '3' => new ChainDemo(),
            '4' => new RaycastDemo(),
            '5' => new TriggerDemo(),
            _ => null
        };

        if (game != null)
        {
            game.Run();
        }

        Console.CursorVisible = true;
        Console.Clear();
        Console.WriteLine("Thanks for trying the Advanced Physics Demo!");
    }
}

abstract class Game
{
    protected PhysicsWorld _world;
    protected AdvancedRenderer _renderer;
    protected bool _running;
    protected const float FixedTimeStep = 1f / 60f;
    protected float _accumulator;
    protected DateTime _lastTime;

    public Game()
    {
        _world = new PhysicsWorld(new Vector2(0, -20f));
        _renderer = new AdvancedRenderer(80, 22, 1.5f); // Smaller viewport to prevent console scrolling
        _running = true;
        _lastTime = DateTime.Now;
    }

    public void Run()
    {
        Initialize();

        while (_running)
        {
            DateTime currentTime = DateTime.Now;
            float deltaTime = (float)(currentTime - _lastTime).TotalSeconds;
            _lastTime = currentTime;

            if (deltaTime > 0.25f)
                deltaTime = 0.25f;

            HandleInput();
            Update(deltaTime);
            Render();

            Thread.Sleep(16);
        }
    }

    protected abstract void Initialize();
    protected abstract void HandleInput();
    protected abstract void CustomRender();

    protected virtual void Update(float deltaTime)
    {
        _accumulator += deltaTime;

        while (_accumulator >= FixedTimeStep)
        {
            _world.Step(FixedTimeStep);
            _accumulator -= FixedTimeStep;
        }
    }

    protected virtual void Render()
    {
        _renderer.Clear();
        CustomRender();
        _renderer.Present();
    }

    protected void CreateGround(float y = 0)
    {
        var ground = new RigidBody(new Vector2(25, y), 0, new BoxShape(100, 2), isStatic: true);
        ground.CollisionLayer = CollisionLayers.Ground;
        _world.AddBody(ground);
    }
}

class RagdollDemo : Game
{
    private RagdollBuilder.Ragdoll? _ragdoll;
    private int _ragdollCount = 0;

    protected override void Initialize()
    {
        CreateGround(-1);

        // Create initial ragdoll
        SpawnRagdoll(new Vector2(15, 15));

        // Create some obstacles
        var box1 = new RigidBody(new Vector2(25, 3), 10, new BoxShape(3, 1));
        _world.AddBody(box1);

        var box2 = new RigidBody(new Vector2(30, 6), 15, new BoxShape(2, 2));
        _world.AddBody(box2);
    }

    private void SpawnRagdoll(Vector2 position)
    {
        _ragdoll = RagdollBuilder.CreateRagdoll(_world, position);
        _ragdollCount++;
    }

    protected override void HandleInput()
    {
        if (!Console.KeyAvailable) return;

        var key = Console.ReadKey(true).Key;

        switch (key)
        {
            case ConsoleKey.Spacebar:
                // Push ragdoll
                _ragdoll?.ApplyImpulseToAll(new Vector2(100, 50));
                break;

            case ConsoleKey.R:
                // Spawn new ragdoll
                SpawnRagdoll(new Vector2(10 + _ragdollCount * 3, 15));
                break;

            case ConsoleKey.Q:
                _running = false;
                break;
        }
    }

    protected override void CustomRender()
    {
        // Draw ground
        _renderer.DrawBox(new Vector2(25, -1), 100, 2, '=', ConsoleColor.DarkGreen);

        // Draw all bodies
        foreach (var body in _world.Bodies)
        {
            if (body.IsStatic) continue;

            ConsoleColor color = body.CollisionLayer == CollisionLayers.Enemy
                ? ConsoleColor.Yellow : ConsoleColor.Gray;

            if (body.Shape is CircleShape circle)
            {
                _renderer.DrawCircle(body.Position, circle.Radius, 'O', color);
            }
            else if (body.Shape is BoxShape box)
            {
                _renderer.DrawBox(body.Position, box.Width, box.Height, '█', color);
            }
        }

        // Draw joints
        foreach (var joint in _world.Joints)
        {
            _renderer.DrawJoint(joint);
        }

        // UI
        _renderer.DrawText(2, 1, "=== RAGDOLL PHYSICS DEMO ===", ConsoleColor.Cyan);
        _renderer.DrawText(2, 2, $"Ragdolls: {_ragdollCount}", ConsoleColor.White);
        _renderer.DrawText(2, 3, $"Bodies: {_world.Bodies.Count}", ConsoleColor.White);
        _renderer.DrawText(2, 4, $"Joints: {_world.Joints.Count}", ConsoleColor.Green);
        _renderer.DrawText(2, 6, "SPACE - Push ragdoll", ConsoleColor.Yellow);
        _renderer.DrawText(2, 7, "R - Spawn new ragdoll", ConsoleColor.Yellow);
        _renderer.DrawText(2, 8, "Q - Quit", ConsoleColor.Yellow);
    }
}

class VehicleDemo : Game
{
    private VehicleBuilder.Vehicle? _vehicle;
    private List<RigidBody> _obstacles = new();

    protected override void Initialize()
    {
        CreateGround(-1);

        // Create terrain with bumps
        for (int i = 0; i < 10; i++)
        {
            float x = 10 + i * 5;
            float height = (i % 2 == 0) ? 1f : 2f;
            var bump = new RigidBody(new Vector2(x, height - 1), 0, new BoxShape(4, height), isStatic: true);
            bump.CollisionLayer = CollisionLayers.Ground;
            _world.AddBody(bump);
        }

        // Create vehicle
        _vehicle = VehicleBuilder.CreateVehicle(_world, new Vector2(10, 5));

        // Create obstacles
        for (int i = 0; i < 5; i++)
        {
            var obstacle = new RigidBody(new Vector2(20 + i * 8, 3), 5, new BoxShape(2, 2));
            _obstacles.Add(obstacle);
            _world.AddBody(obstacle);
        }
    }

    protected override void HandleInput()
    {
        if (!Console.KeyAvailable) return;

        var key = Console.ReadKey(true).Key;

        switch (key)
        {
            case ConsoleKey.LeftArrow:
                _vehicle?.ApplyAcceleration(-500f);
                break;

            case ConsoleKey.RightArrow:
                _vehicle?.ApplyAcceleration(500f);
                break;

            case ConsoleKey.UpArrow:
                _vehicle?.Body.ApplyImpulse(new Vector2(0, 100));
                break;

            case ConsoleKey.Q:
                _running = false;
                break;
        }
    }

    protected override void CustomRender()
    {
        // Draw terrain
        foreach (var body in _world.Bodies)
        {
            if (body.CollisionLayer == CollisionLayers.Ground && body.Shape is BoxShape box)
            {
                _renderer.DrawBox(body.Position, box.Width, box.Height, '▓', ConsoleColor.DarkGreen);
            }
        }

        // Draw obstacles
        foreach (var obstacle in _obstacles)
        {
            if (obstacle.Shape is BoxShape box)
            {
                _renderer.DrawBox(obstacle.Position, box.Width, box.Height, '■', ConsoleColor.Red);
            }
        }

        // Draw vehicle
        if (_vehicle != null)
        {
            if (_vehicle.Body.Shape is BoxShape bodyBox)
            {
                _renderer.DrawBox(_vehicle.Body.Position, bodyBox.Width, bodyBox.Height, '█', ConsoleColor.Cyan);
            }

            if (_vehicle.FrontWheel.Shape is CircleShape wheel1)
            {
                _renderer.DrawCircle(_vehicle.FrontWheel.Position, wheel1.Radius, 'O', ConsoleColor.White);
            }

            if (_vehicle.BackWheel.Shape is CircleShape wheel2)
            {
                _renderer.DrawCircle(_vehicle.BackWheel.Position, wheel2.Radius, 'O', ConsoleColor.White);
            }

            // Draw suspension springs
            foreach (var joint in _vehicle.Joints)
            {
                _renderer.DrawJoint(joint);
            }
        }

        // UI
        _renderer.DrawText(2, 1, "=== VEHICLE DEMO ===", ConsoleColor.Cyan);
        _renderer.DrawText(2, 2, "Spring suspension with revolute joints", ConsoleColor.Gray);
        _renderer.DrawText(2, 4, "LEFT/RIGHT - Accelerate", ConsoleColor.Yellow);
        _renderer.DrawText(2, 5, "UP - Jump", ConsoleColor.Yellow);
        _renderer.DrawText(2, 6, "Q - Quit", ConsoleColor.Yellow);

        if (_vehicle != null)
        {
            _renderer.DrawText(2, 8, $"Velocity: {_vehicle.Body.Velocity.Length:F1}", ConsoleColor.White);
        }
    }
}

class ChainDemo : Game
{
    private List<RigidBody> _chain = new();
    private List<RigidBody> _bridge = new();

    protected override void Initialize()
    {
        CreateGround(-1);

        // Create hanging chain
        _chain = ChainBuilder.CreateChain(_world, new Vector2(10, 15), 12, 1.0f, isRope: false);

        // Create rope
        ChainBuilder.CreateChain(_world, new Vector2(20, 18), 15, 0.8f, isRope: true);

        // Create bridge
        _bridge = ChainBuilder.CreateBridge(_world, new Vector2(25, 10), 8, 2.0f);

        // Add weight to bridge
        var weight = new RigidBody(new Vector2(33, 12), 20, new CircleShape(1.0f));
        weight.Friction = 0.8f;
        _world.AddBody(weight);
    }

    protected override void HandleInput()
    {
        if (!Console.KeyAvailable) return;

        var key = Console.ReadKey(true).Key;

        switch (key)
        {
            case ConsoleKey.Spacebar:
                // Add sphere to interact with chains
                var sphere = new RigidBody(new Vector2(15, 20), 5, new CircleShape(0.8f));
                sphere.Restitution = 0.7f;
                _world.AddBody(sphere);
                break;

            case ConsoleKey.Q:
                _running = false;
                break;
        }
    }

    protected override void CustomRender()
    {
        // Draw ground
        _renderer.DrawBox(new Vector2(25, -1), 100, 2, '=', ConsoleColor.DarkGreen);

        // Draw all bodies
        foreach (var body in _world.Bodies)
        {
            if (body.IsStatic && body.CollisionLayer != CollisionLayers.Ground) continue;

            ConsoleColor color = ConsoleColor.Gray;

            if (body.Shape is CircleShape circle)
            {
                _renderer.DrawCircle(body.Position, circle.Radius, 'O', color);
            }
            else if (body.Shape is BoxShape box)
            {
                _renderer.DrawBox(body.Position, box.Width, box.Height, '▓', color);
            }
        }

        // Draw all joints (chains, ropes, bridge)
        foreach (var joint in _world.Joints)
        {
            _renderer.DrawJoint(joint);
        }

        // UI
        _renderer.DrawText(2, 1, "=== CHAINS & BRIDGE DEMO ===", ConsoleColor.Cyan);
        _renderer.DrawText(2, 2, "Distance joints, rope joints", ConsoleColor.Gray);
        _renderer.DrawText(2, 4, $"Joints: {_world.Joints.Count}", ConsoleColor.Green);
        _renderer.DrawText(2, 5, $"Bodies: {_world.Bodies.Count}", ConsoleColor.White);
        _renderer.DrawText(2, 7, "SPACE - Drop sphere", ConsoleColor.Yellow);
        _renderer.DrawText(2, 8, "Q - Quit", ConsoleColor.Yellow);
    }
}

class RaycastDemo : Game
{
    private Vector2 _rayOrigin = new Vector2(5, 10);
    private Vector2 _rayDirection = new Vector2(1, 0);
    private RaycastHit _lastHit;
    private List<RigidBody> _targets = new();

    protected override void Initialize()
    {
        CreateGround(-1);

        // Create targets
        for (int i = 0; i < 8; i++)
        {
            float x = 15 + (i % 4) * 7;
            float y = 3 + (i / 4) * 8;

            Shape shape = i % 2 == 0 ? (Shape)new CircleShape(1.5f) : new BoxShape(2, 2);
            var target = new RigidBody(new Vector2(x, y), 10, shape);
            target.CollisionLayer = CollisionLayers.Enemy;
            target.Restitution = 0.5f;

            _targets.Add(target);
            _world.AddBody(target);
        }
    }

    protected override void HandleInput()
    {
        if (!Console.KeyAvailable) return;

        var key = Console.ReadKey(true).Key;

        switch (key)
        {
            case ConsoleKey.W:
                _rayDirection = new Vector2(_rayDirection.X, _rayDirection.Y + 0.1f).Normalized;
                break;

            case ConsoleKey.S:
                _rayDirection = new Vector2(_rayDirection.X, _rayDirection.Y - 0.1f).Normalized;
                break;

            case ConsoleKey.Spacebar:
                // Shoot - apply impulse to hit object
                var ray = new Ray(_rayOrigin, _rayDirection, 50f);
                var hit = _world.Raycast(ray, CollisionLayers.Enemy);

                if (hit.Hit)
                {
                    hit.Body.ApplyImpulseAtPoint(_rayDirection * 50f, hit.Point);
                }
                break;

            case ConsoleKey.Q:
                _running = false;
                break;
        }
    }

    protected override void Update(float deltaTime)
    {
        base.Update(deltaTime);

        // Update raycast
        var ray = new Ray(_rayOrigin, _rayDirection, 50f);
        _lastHit = _world.Raycast(ray, CollisionLayers.Enemy);
    }

    protected override void CustomRender()
    {
        // Draw ground
        _renderer.DrawBox(new Vector2(25, -1), 100, 2, '=', ConsoleColor.DarkGreen);

        // Draw targets
        foreach (var target in _targets)
        {
            if (target.Shape is CircleShape circle)
            {
                _renderer.DrawCircle(target.Position, circle.Radius, 'O', ConsoleColor.Red);
            }
            else if (target.Shape is BoxShape box)
            {
                _renderer.DrawBox(target.Position, box.Width, box.Height, '■', ConsoleColor.Red);
            }
        }

        // Draw ray
        var ray = new Ray(_rayOrigin, _rayDirection, 50f);
        _renderer.DrawRaycast(ray, _lastHit);

        // Draw ray origin
        _renderer.DrawCircle(_rayOrigin, 0.5f, '@', ConsoleColor.Cyan);

        // UI
        _renderer.DrawText(2, 1, "=== RAYCASTING DEMO ===", ConsoleColor.Cyan);
        _renderer.DrawText(2, 2, "Raycast physics queries", ConsoleColor.Gray);
        _renderer.DrawText(2, 4, "W/S - Aim", ConsoleColor.Yellow);
        _renderer.DrawText(2, 5, "SPACE - Shoot", ConsoleColor.Yellow);
        _renderer.DrawText(2, 6, "Q - Quit", ConsoleColor.Yellow);

        if (_lastHit.Hit)
        {
            _renderer.DrawText(2, 8, $"HIT! Distance: {_lastHit.Distance:F2}", ConsoleColor.Green);
        }
        else
        {
            _renderer.DrawText(2, 8, "No hit", ConsoleColor.DarkGray);
        }

        _renderer.DrawText(2, 9, $"Angle: {MathF.Atan2(_rayDirection.Y, _rayDirection.X) * 180 / MathF.PI:F0}°", ConsoleColor.White);
    }
}

class TriggerDemo : Game
{
    private List<RigidBody> _triggers = new();
    private List<RigidBody> _balls = new();
    private int _triggerEnterCount = 0;
    private string _lastEvent = "";

    protected override void Initialize()
    {
        CreateGround(-1);

        // Create trigger zones
        var trigger1 = new RigidBody(new Vector2(15, 5), 0, new CircleShape(3f), isStatic: true);
        trigger1.IsTrigger = true;
        trigger1.CollisionLayer = CollisionLayers.Trigger;
        _triggers.Add(trigger1);
        _world.AddBody(trigger1);

        var trigger2 = new RigidBody(new Vector2(30, 8), 0, new BoxShape(5, 5), isStatic: true);
        trigger2.IsTrigger = true;
        trigger2.CollisionLayer = CollisionLayers.Trigger;
        _triggers.Add(trigger2);
        _world.AddBody(trigger2);

        // Subscribe to events
        _world.OnTriggerEnter += OnTriggerEnter;
        _world.OnTriggerExit += OnTriggerExit;
        _world.OnCollisionEnter += OnCollisionEnter;

        // Create walls
        var wall1 = new RigidBody(new Vector2(0, 10), 0, new BoxShape(1, 30), isStatic: true);
        _world.AddBody(wall1);

        var wall2 = new RigidBody(new Vector2(50, 10), 0, new BoxShape(1, 30), isStatic: true);
        _world.AddBody(wall2);
    }

    private void OnTriggerEnter(object? sender, CollisionEventArgs e)
    {
        _triggerEnterCount++;
        _lastEvent = $"Trigger ENTER #{_triggerEnterCount}";
    }

    private void OnTriggerExit(object? sender, CollisionEventArgs e)
    {
        _lastEvent = "Trigger EXIT";
    }

    private void OnCollisionEnter(object? sender, CollisionEventArgs e)
    {
        if (!e.BodyA.IsStatic && !e.BodyB.IsStatic)
        {
            _lastEvent = "Ball collision!";
        }
    }

    protected override void HandleInput()
    {
        if (!Console.KeyAvailable) return;

        var key = Console.ReadKey(true).Key;

        switch (key)
        {
            case ConsoleKey.Spacebar:
                // Spawn ball
                var ball = new RigidBody(new Vector2(5 + _balls.Count * 2, 15), 3, new CircleShape(0.8f));
                ball.Restitution = 0.6f;
                ball.CollisionLayer = CollisionLayers.Default;
                ball.CollisionMask = CollisionLayers.Everything;
                _balls.Add(ball);
                _world.AddBody(ball);
                break;

            case ConsoleKey.Q:
                _running = false;
                break;
        }
    }

    protected override void CustomRender()
    {
        // Draw ground
        _renderer.DrawBox(new Vector2(25, -1), 100, 2, '=', ConsoleColor.DarkGreen);

        // Draw trigger zones (outlined)
        foreach (var trigger in _triggers)
        {
            if (trigger.Shape is CircleShape circle)
            {
                _renderer.DrawCircle(trigger.Position, circle.Radius, '░', ConsoleColor.DarkCyan, filled: false);
            }
            else if (trigger.Shape is BoxShape box)
            {
                _renderer.DrawBox(trigger.Position, box.Width, box.Height, '░', ConsoleColor.DarkCyan, filled: false);
            }
        }

        // Draw balls
        foreach (var ball in _balls)
        {
            if (ball.Shape is CircleShape circle)
            {
                _renderer.DrawCircle(ball.Position, circle.Radius, 'O', ConsoleColor.Yellow);
            }
        }

        // Draw walls
        foreach (var body in _world.Bodies)
        {
            if (body.IsStatic && !body.IsTrigger && body.CollisionLayer != CollisionLayers.Ground)
            {
                if (body.Shape is BoxShape box)
                {
                    _renderer.DrawBox(body.Position, box.Width, box.Height, '█', ConsoleColor.DarkGray);
                }
            }
        }

        // UI
        _renderer.DrawText(2, 1, "=== TRIGGER ZONES & EVENTS ===", ConsoleColor.Cyan);
        _renderer.DrawText(2, 2, "Collision events, trigger zones", ConsoleColor.Gray);
        _renderer.DrawText(2, 4, "SPACE - Spawn ball", ConsoleColor.Yellow);
        _renderer.DrawText(2, 5, "Q - Quit", ConsoleColor.Yellow);
        _renderer.DrawText(2, 7, $"Balls: {_balls.Count}", ConsoleColor.White);
        _renderer.DrawText(2, 8, $"Trigger enters: {_triggerEnterCount}", ConsoleColor.Green);
        _renderer.DrawText(2, 9, $"Last event: {_lastEvent}", ConsoleColor.Magenta);
    }
}
