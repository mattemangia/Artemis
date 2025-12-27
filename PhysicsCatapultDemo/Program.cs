using ArtemisEngine;
using PhysicsCatapultDemo;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== Physics Catapult Demo ===");
        Console.WriteLine("A demonstration of the Artemis Physics Engine");
        Console.WriteLine();
        Console.WriteLine("Controls:");
        Console.WriteLine("  W/S - Adjust launch angle up/down");
        Console.WriteLine("  A/D - Adjust launch power down/up");
        Console.WriteLine("  SPACE - Launch projectile");
        Console.WriteLine("  R - Reset level");
        Console.WriteLine("  1/2/3 - Load level 1/2/3");
        Console.WriteLine("  Q - Quit");
        Console.WriteLine();
        Console.WriteLine("Press any key to start...");
        Console.ReadKey(true);

        Console.CursorVisible = false;
        Console.Clear();

        var game = new Game();
        game.Run();

        Console.CursorVisible = true;
        Console.Clear();
        Console.WriteLine("Thanks for playing!");
    }
}

class Game
{
    private PhysicsWorld _world;
    private Catapult _catapult;
    private List<GameObject> _gameObjects;
    private Renderer _renderer;
    private bool _running;
    private int _currentLevel;
    private const float FixedTimeStep = 1f / 60f;
    private float _accumulator;
    private DateTime _lastTime;

    public Game()
    {
        _world = new PhysicsWorld(new Vector2(0, -20f)); // Gravity
        _catapult = new Catapult(new Vector2(5, 5));
        _gameObjects = new List<GameObject>();
        _renderer = new Renderer(80, 22, 1.5f); // Smaller viewport to prevent console scrolling
        _running = true;
        _currentLevel = 1;
        _lastTime = DateTime.Now;

        LoadLevel(1);
    }

    public void Run()
    {
        while (_running)
        {
            DateTime currentTime = DateTime.Now;
            float deltaTime = (float)(currentTime - _lastTime).TotalSeconds;
            _lastTime = currentTime;

            // Cap delta time to prevent spiral of death
            if (deltaTime > 0.25f)
                deltaTime = 0.25f;

            HandleInput();
            Update(deltaTime);
            Render();

            Thread.Sleep(16); // ~60 FPS
        }
    }

    private void HandleInput()
    {
        if (!Console.KeyAvailable)
            return;

        var key = Console.ReadKey(true).Key;

        switch (key)
        {
            case ConsoleKey.W:
                _catapult.SetAimAngle(GetAimAngle() + 5);
                break;

            case ConsoleKey.S:
                _catapult.SetAimAngle(GetAimAngle() - 5);
                break;

            case ConsoleKey.A:
                _catapult.AdjustPower(-50);
                break;

            case ConsoleKey.D:
                _catapult.AdjustPower(50);
                break;

            case ConsoleKey.Spacebar:
                LaunchProjectile();
                break;

            case ConsoleKey.R:
                LoadLevel(_currentLevel);
                break;

            case ConsoleKey.D1:
                LoadLevel(1);
                break;

            case ConsoleKey.D2:
                LoadLevel(2);
                break;

            case ConsoleKey.D3:
                LoadLevel(3);
                break;

            case ConsoleKey.Q:
                _running = false;
                break;
        }
    }

    private float GetAimAngle()
    {
        return MathF.Atan2(-_catapult.AimDirection.Y, _catapult.AimDirection.X) * 180f / MathF.PI;
    }

    private void LaunchProjectile()
    {
        var projectile = _catapult.LaunchProjectile(0.6f, 2.0f);
        _gameObjects.Add(projectile);
        _world.AddBody(projectile.Body);
    }

    private void LoadLevel(int level)
    {
        // Clear existing objects
        foreach (var obj in _gameObjects)
        {
            _world.RemoveBody(obj.Body);
        }
        _gameObjects.Clear();

        _currentLevel = level;

        // Load level
        switch (level)
        {
            case 1:
                _gameObjects = LevelBuilder.BuildLevel1();
                break;
            case 2:
                _gameObjects = LevelBuilder.BuildLevel2();
                break;
            case 3:
                _gameObjects = LevelBuilder.BuildLevel3();
                break;
            default:
                _gameObjects = LevelBuilder.BuildLevel1();
                break;
        }

        // Add all bodies to physics world
        foreach (var obj in _gameObjects)
        {
            _world.AddBody(obj.Body);
        }
    }

    private void Update(float deltaTime)
    {
        _accumulator += deltaTime;

        while (_accumulator >= FixedTimeStep)
        {
            _world.Step(FixedTimeStep);
            ProcessCollisionDamage();
            RemoveDestroyedObjects();
            _accumulator -= FixedTimeStep;
        }
    }

    private void ProcessCollisionDamage()
    {
        // Simple damage system based on velocity
        foreach (var obj in _gameObjects)
        {
            if (obj.Type == GameObjectType.Block && obj.Body.Velocity.Length > 5f)
            {
                float damage = obj.Body.Velocity.Length * 2f;
                obj.TakeDamage(damage * FixedTimeStep);
            }
        }
    }

    private void RemoveDestroyedObjects()
    {
        var toRemove = _gameObjects.Where(obj => obj.IsDestroyed).ToList();

        foreach (var obj in toRemove)
        {
            _world.RemoveBody(obj.Body);
            _gameObjects.Remove(obj);
        }

        // Also remove projectiles that are out of bounds or too slow
        var projectilesToRemove = _gameObjects
            .Where(obj => obj.Type == GameObjectType.Projectile &&
                         (obj.Body.Position.Y < -10 ||
                          obj.Body.Position.X > 60 ||
                          obj.Body.Position.X < -10 ||
                          (obj.Body.Velocity.Length < 0.5f && obj.Body.Position.Y < 15)))
            .ToList();

        foreach (var proj in projectilesToRemove)
        {
            _world.RemoveBody(proj.Body);
            _gameObjects.Remove(proj);
        }
    }

    private void Render()
    {
        _renderer.Clear();

        // Draw ground
        foreach (var obj in _gameObjects.Where(o => o.Type == GameObjectType.Ground))
        {
            if (obj.Body.Shape is BoxShape box)
            {
                _renderer.DrawBox(obj.Body.Position, box.Width, box.Height, '=', ConsoleColor.DarkGreen);
            }
        }

        // Draw blocks
        foreach (var obj in _gameObjects.Where(o => o.Type == GameObjectType.Block))
        {
            if (obj.Body.Shape is BoxShape box && obj.Material != null)
            {
                char c = obj.Material.Name switch
                {
                    "Wood" => '#',
                    "Stone" => '█',
                    "Glass" => '▒',
                    "Metal" => '▓',
                    _ => '■'
                };

                // Adjust color based on health
                ConsoleColor color = obj.Material.Color;
                if (obj.Health < obj.MaxHealth * 0.3f)
                {
                    color = ConsoleColor.Red;
                }
                else if (obj.Health < obj.MaxHealth * 0.6f)
                {
                    color = ConsoleColor.DarkYellow;
                }

                _renderer.DrawBox(obj.Body.Position, box.Width, box.Height, c, color);
            }
        }

        // Draw projectiles
        foreach (var obj in _gameObjects.Where(o => o.Type == GameObjectType.Projectile))
        {
            if (obj.Body.Shape is CircleShape circle)
            {
                _renderer.DrawCircle(obj.Body.Position, circle.Radius, 'O', ConsoleColor.Red);
            }
        }

        // Draw catapult
        _renderer.DrawCircle(_catapult.Position, 0.8f, 'C', ConsoleColor.DarkMagenta);
        _renderer.DrawLine(_catapult.Position, _catapult.GetAimEndPoint(), '-', ConsoleColor.Magenta);

        // Draw UI
        _renderer.DrawText(2, 1, $"Level: {_currentLevel}", ConsoleColor.White);
        _renderer.DrawText(2, 2, $"Angle: {GetAimAngle():F0}°", ConsoleColor.White);
        _renderer.DrawText(2, 3, $"Power: {_catapult.Power:F0}", ConsoleColor.White);
        _renderer.DrawText(2, 4, $"Blocks: {_gameObjects.Count(o => o.Type == GameObjectType.Block)}", ConsoleColor.White);

        int projectileCount = _gameObjects.Count(o => o.Type == GameObjectType.Projectile);
        if (projectileCount > 0)
        {
            _renderer.DrawText(2, 5, $"Projectiles: {projectileCount}", ConsoleColor.Yellow);
        }

        // Check win condition
        if (_gameObjects.Count(o => o.Type == GameObjectType.Block) == 0)
        {
            _renderer.DrawText(40, 10, "*** LEVEL COMPLETE ***", ConsoleColor.Green);
            _renderer.DrawText(40, 11, "Press 1/2/3 for next level", ConsoleColor.Green);
        }

        _renderer.Present();
    }
}
