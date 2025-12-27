using Artemis.Graphics;
using ArtemisEngine;
using OpenTK.Mathematics;
using OpenTK.Windowing.GraphicsLibraryFramework;
using PhysicsCatapultDemo;

namespace PhysicsCatapultDemo.OpenTK;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== Physics Catapult Demo - OpenTK Edition ===");
        Console.WriteLine("Starting graphical demo...");

        using var window = new CatapultWindow(1280, 720, "Artemis Physics Catapult");
        window.Run();
    }
}

public class CatapultWindow : GraphicsWindow
{
    private Catapult _catapult = null!;
    private List<GameObject> _gameObjects = new();
    private int _currentLevel = 1;

    // Material colors
    private static readonly Dictionary<string, Color4> MaterialColors = new()
    {
        { "Wood", new Color4(0.65f, 0.45f, 0.25f, 1f) },
        { "Stone", new Color4(0.5f, 0.5f, 0.55f, 1f) },
        { "Glass", new Color4(0.6f, 0.8f, 0.9f, 0.7f) },
        { "Metal", new Color4(0.7f, 0.7f, 0.75f, 1f) },
    };

    public CatapultWindow(int width, int height, string title)
        : base(width, height, title)
    {
    }

    protected override void Initialize()
    {
        World = new PhysicsWorld(new ArtemisEngine.Vector2(0, -20f));
        _catapult = new Catapult(new ArtemisEngine.Vector2(5, 8));

        // Camera position
        CameraPosition = new Vector2d(30, 15);
        CameraZoom = 25;

        LoadLevel(1);
    }

    private void LoadLevel(int level)
    {
        // Clear existing objects
        foreach (var obj in _gameObjects)
        {
            World.RemoveBody(obj.Body);
        }
        _gameObjects.Clear();

        _currentLevel = level;

        // Load level
        _gameObjects = level switch
        {
            1 => LevelBuilder.BuildLevel1(),
            2 => LevelBuilder.BuildLevel2(),
            3 => LevelBuilder.BuildLevel3(),
            _ => LevelBuilder.BuildLevel1()
        };

        // Add all bodies to physics world
        foreach (var obj in _gameObjects)
        {
            World.AddBody(obj.Body);
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
        World.AddBody(projectile.Body);
    }

    protected override void HandleInput()
    {
        // Aim controls
        if (KeyboardState.IsKeyDown(Keys.W))
            _catapult.SetAimAngle(GetAimAngle() + 1);
        if (KeyboardState.IsKeyDown(Keys.S))
            _catapult.SetAimAngle(GetAimAngle() - 1);

        // Power controls
        if (KeyboardState.IsKeyDown(Keys.A))
            _catapult.AdjustPower(-5);
        if (KeyboardState.IsKeyDown(Keys.D))
            _catapult.AdjustPower(5);

        // Launch
        if (KeyboardState.IsKeyPressed(Keys.Space))
            LaunchProjectile();

        // Level selection
        if (KeyboardState.IsKeyPressed(Keys.D1))
            LoadLevel(1);
        if (KeyboardState.IsKeyPressed(Keys.D2))
            LoadLevel(2);
        if (KeyboardState.IsKeyPressed(Keys.D3))
            LoadLevel(3);

        // Reset level
        if (KeyboardState.IsKeyPressed(Keys.R))
            LoadLevel(_currentLevel);
    }

    protected override void UpdatePhysics(float deltaTime)
    {
        // Process collision damage
        foreach (var obj in _gameObjects)
        {
            if (obj.Type == GameObjectType.Block && obj.Body.Velocity.Length > 5f)
            {
                float damage = obj.Body.Velocity.Length * 2f;
                obj.TakeDamage(damage * deltaTime);
            }
        }

        // Remove destroyed objects
        var toRemove = _gameObjects.Where(obj => obj.IsDestroyed).ToList();
        foreach (var obj in toRemove)
        {
            World.RemoveBody(obj.Body);
            _gameObjects.Remove(obj);
        }

        // Remove out-of-bounds projectiles
        var projectilesToRemove = _gameObjects
            .Where(obj => obj.Type == GameObjectType.Projectile &&
                         (obj.Body.Position.Y < -10 ||
                          obj.Body.Position.X > 80 ||
                          obj.Body.Position.X < -10 ||
                          (obj.Body.Velocity.Length < 0.5f && obj.Body.Position.Y < 15)))
            .ToList();

        foreach (var proj in projectilesToRemove)
        {
            World.RemoveBody(proj.Body);
            _gameObjects.Remove(proj);
        }
    }

    protected override void Render()
    {
        // Draw ground
        foreach (var obj in _gameObjects.Where(o => o.Type == GameObjectType.Ground))
        {
            if (obj.Body.Shape is BoxShape box)
            {
                DrawBox(obj.Body.Position, box.Width, box.Height, obj.Body.Rotation,
                    new Color4(0.2f, 0.5f, 0.2f, 1f));
            }
        }

        // Draw blocks
        foreach (var obj in _gameObjects.Where(o => o.Type == GameObjectType.Block))
        {
            if (obj.Body.Shape is BoxShape box && obj.Material != null)
            {
                Color4 color = MaterialColors.GetValueOrDefault(obj.Material.Name, Color4.Gray);

                // Adjust color based on health
                float healthRatio = obj.Health / obj.MaxHealth;
                if (healthRatio < 0.3f)
                {
                    color = new Color4(0.9f, 0.2f, 0.2f, color.A);
                }
                else if (healthRatio < 0.6f)
                {
                    color = new Color4(0.9f, 0.7f, 0.2f, color.A);
                }

                DrawBox(obj.Body.Position, box.Width, box.Height, obj.Body.Rotation, color);

                // Draw outline
                DrawBox(obj.Body.Position, box.Width, box.Height, obj.Body.Rotation,
                    new Color4(0f, 0f, 0f, 0.5f), false);
            }
        }

        // Draw projectiles
        foreach (var obj in _gameObjects.Where(o => o.Type == GameObjectType.Projectile))
        {
            if (obj.Body.Shape is CircleShape circle)
            {
                // Glow effect
                DrawCircle(obj.Body.Position, circle.Radius * 1.5f,
                    new Color4(1f, 0.3f, 0.1f, 0.3f), true);
                // Main projectile
                DrawCircle(obj.Body.Position, circle.Radius,
                    new Color4(0.9f, 0.2f, 0.1f, 1f), true);
            }
        }

        // Draw catapult
        DrawCircle(_catapult.Position, 1.2f, new Color4(0.5f, 0.2f, 0.4f, 1f), true);
        DrawCircle(_catapult.Position, 1.2f, new Color4(0.3f, 0.1f, 0.3f, 1f), false);

        // Draw aim line
        ArtemisEngine.Vector2 aimEnd = _catapult.Position + _catapult.AimDirection * 8f;
        DrawLine(_catapult.Position, aimEnd, new Color4(1f, 1f, 0.2f, 0.8f), 2f);

        // Draw power indicator
        float powerRatio = (_catapult.Power - _catapult.MinPower) / (_catapult.MaxPower - _catapult.MinPower);
        ArtemisEngine.Vector2 powerEnd = _catapult.Position + _catapult.AimDirection * (3f + powerRatio * 5f);
        DrawLine(_catapult.Position, powerEnd, new Color4(1f, powerRatio, 0f, 1f), 3f);

        // Draw trajectory preview (simple arc)
        DrawTrajectoryPreview();

        // Check win condition
        int blockCount = _gameObjects.Count(o => o.Type == GameObjectType.Block);
        if (blockCount == 0)
        {
            // Draw level complete indicator
            DrawCircle(new ArtemisEngine.Vector2((float)CameraPosition.X, (float)CameraPosition.Y + 5f), 3f,
                new Color4(0.2f, 0.9f, 0.2f, 0.5f), true);
        }
    }

    private void DrawTrajectoryPreview()
    {
        const int points = 30;
        const float timeStep = 0.05f;

        ArtemisEngine.Vector2 pos = _catapult.Position + _catapult.AimDirection * 2f;
        ArtemisEngine.Vector2 vel = _catapult.AimDirection * _catapult.Power;
        ArtemisEngine.Vector2 gravity = new ArtemisEngine.Vector2(0, -20f);

        List<ArtemisEngine.Vector2> trajectory = new() { pos };

        for (int i = 0; i < points; i++)
        {
            vel = vel + gravity * timeStep;
            pos = pos + vel * timeStep;

            if (pos.Y < 0) break;

            trajectory.Add(pos);
        }

        if (trajectory.Count > 1)
        {
            DrawTrail(trajectory,
                new Color4(1f, 1f, 0.5f, 0.1f),
                new Color4(1f, 0.5f, 0.2f, 0.5f));
        }
    }

    protected override void DrawUI()
    {
        int blockCount = _gameObjects.Count(o => o.Type == GameObjectType.Block);
        int projectileCount = _gameObjects.Count(o => o.Type == GameObjectType.Projectile);

        string status = blockCount == 0 ? "LEVEL COMPLETE! " : "";

        Title = $"Artemis Catapult | Level: {_currentLevel} | " +
                $"Angle: {GetAimAngle():F0}° | " +
                $"Power: {_catapult.Power:F0} | " +
                $"Blocks: {blockCount} | " +
                $"Projectiles: {projectileCount} | " +
                status +
                (IsPaused ? "[PAUSED] " : "") +
                "[SPACE]Fire [W/S]Aim [A/D]Power [1-3]Levels [R]Reset [ESC]Quit";
    }
}
