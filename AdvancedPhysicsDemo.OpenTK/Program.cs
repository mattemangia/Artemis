using Artemis.Graphics;
using ArtemisEngine;
using OpenTK.Mathematics;
using OpenTK.Windowing.GraphicsLibraryFramework;
using AdvancedPhysicsDemo;

namespace AdvancedPhysicsDemo.OpenTK;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== Advanced Physics Demo - OpenTK Edition ===");
        Console.WriteLine("Demonstrating joints, ragdolls, vehicles, and chains");

        using var window = new AdvancedPhysicsWindow(1280, 720, "Artemis Advanced Physics Demo");
        window.Run();
    }
}

public class AdvancedPhysicsWindow : GraphicsWindow
{
    private RagdollBuilder.Ragdoll? _ragdoll;
    private VehicleBuilder.Vehicle? _vehicle;
    private List<RigidBody> _chain = new();
    private List<RigidBody> _bridge = new();
    private RigidBody? _ground;

    // Demo mode
    private int _currentDemo = 0;
    private const int DemoCount = 4;
    private string[] _demoNames = { "Ragdoll", "Vehicle", "Chain", "Bridge" };

    public AdvancedPhysicsWindow(int width, int height, string title)
        : base(width, height, title)
    {
    }

    protected override void Initialize()
    {
        World = new PhysicsWorld(new ArtemisEngine.Vector2(0, -20f));
        World.UseAdvancedSolver = true;

        // Camera settings
        CameraPosition = new Vector2d(15, 10);
        CameraZoom = 25;

        // Create ground
        CreateGround();

        // Load first demo
        LoadDemo(0);
    }

    private void CreateGround()
    {
        _ground = new RigidBody(new ArtemisEngine.Vector2(20, 0), 0, new BoxShape(60, 1), isStatic: true);
        _ground.Friction = 0.8f;
        _ground.CollisionLayer = CollisionLayers.Static;
        World.AddBody(_ground);
    }

    private void ClearDemo()
    {
        // Remove all non-ground bodies
        var bodiesToRemove = World.Bodies.Where(b => b != _ground).ToList();
        foreach (var body in bodiesToRemove)
        {
            World.RemoveBody(body);
        }

        // Clear joints
        foreach (var joint in World.Joints.ToList())
        {
            World.RemoveJoint(joint);
        }

        _ragdoll = null;
        _vehicle = null;
        _chain.Clear();
        _bridge.Clear();
    }

    private void LoadDemo(int demoIndex)
    {
        ClearDemo();
        _currentDemo = demoIndex;

        switch (demoIndex)
        {
            case 0:
                LoadRagdollDemo();
                break;
            case 1:
                LoadVehicleDemo();
                break;
            case 2:
                LoadChainDemo();
                break;
            case 3:
                LoadBridgeDemo();
                break;
        }

        CameraPosition = new Vector2d(15, 10);
    }

    private void LoadRagdollDemo()
    {
        _ragdoll = RagdollBuilder.CreateRagdoll(World, new ArtemisEngine.Vector2(10, 5));

        // Add some platforms
        var platform1 = new RigidBody(new ArtemisEngine.Vector2(5, 3), 0, new BoxShape(6, 0.5f), isStatic: true);
        platform1.CollisionLayer = CollisionLayers.Static;
        World.AddBody(platform1);

        var platform2 = new RigidBody(new ArtemisEngine.Vector2(20, 6), 0, new BoxShape(6, 0.5f), isStatic: true);
        platform2.Rotation = -0.2f;
        platform2.CollisionLayer = CollisionLayers.Static;
        World.AddBody(platform2);
    }

    private void LoadVehicleDemo()
    {
        _vehicle = VehicleBuilder.CreateVehicle(World, new ArtemisEngine.Vector2(5, 5));

        // Create ramp
        var ramp = new RigidBody(new ArtemisEngine.Vector2(25, 2), 0, new BoxShape(10, 0.5f), isStatic: true);
        ramp.Rotation = 0.3f;
        ramp.CollisionLayer = CollisionLayers.Static;
        World.AddBody(ramp);

        // Create obstacles
        for (int i = 0; i < 5; i++)
        {
            var box = new RigidBody(new ArtemisEngine.Vector2(35 + i * 1.5f, 1 + i * 0.5f), 2f, new BoxShape(1, 1));
            box.Friction = 0.5f;
            World.AddBody(box);
        }
    }

    private void LoadChainDemo()
    {
        _chain = ChainBuilder.CreateChain(World, new ArtemisEngine.Vector2(10, 15), 10, 1.2f, isRope: false);

        // Add a swinging weight
        var weight = new RigidBody(new ArtemisEngine.Vector2(10, 3), 10f, new CircleShape(1.5f));
        weight.Friction = 0.6f;
        World.AddBody(weight);

        var lastChainBody = _chain.Last();
        var ropeJoint = new RopeJoint(lastChainBody, weight, 2f);
        World.AddJoint(ropeJoint);
    }

    private void LoadBridgeDemo()
    {
        _bridge = ChainBuilder.CreateBridge(World, new ArtemisEngine.Vector2(5, 8), 8, 2f);

        // Add some balls to bounce on the bridge
        for (int i = 0; i < 3; i++)
        {
            var ball = new RigidBody(new ArtemisEngine.Vector2(10 + i * 3, 12 + i * 2), 3f, new CircleShape(0.8f));
            ball.Friction = 0.4f;
            ball.Restitution = 0.5f;
            World.AddBody(ball);
        }
    }

    protected override void HandleInput()
    {
        // Demo switching
        if (KeyboardState.IsKeyPressed(Keys.D1)) LoadDemo(0);
        if (KeyboardState.IsKeyPressed(Keys.D2)) LoadDemo(1);
        if (KeyboardState.IsKeyPressed(Keys.D3)) LoadDemo(2);
        if (KeyboardState.IsKeyPressed(Keys.D4)) LoadDemo(3);

        // Reload current demo
        if (KeyboardState.IsKeyPressed(Keys.R))
            LoadDemo(_currentDemo);

        // Demo-specific controls
        switch (_currentDemo)
        {
            case 0: // Ragdoll
                if (KeyboardState.IsKeyPressed(Keys.Space) && _ragdoll != null)
                {
                    _ragdoll.ApplyImpulseToAll(new ArtemisEngine.Vector2(50, 100));
                }
                break;

            case 1: // Vehicle
                if (_vehicle != null)
                {
                    if (KeyboardState.IsKeyDown(Keys.A))
                        _vehicle.ApplyAcceleration(-50);
                    if (KeyboardState.IsKeyDown(Keys.D))
                        _vehicle.ApplyAcceleration(50);
                    if (KeyboardState.IsKeyPressed(Keys.Space))
                        _vehicle.Body.ApplyImpulse(new ArtemisEngine.Vector2(0, 200));
                }
                break;

            case 2: // Chain
                if (KeyboardState.IsKeyPressed(Keys.Space) && _chain.Count > 1)
                {
                    var lastBody = _chain.Last();
                    lastBody.ApplyImpulse(new ArtemisEngine.Vector2(100, 0));
                }
                break;

            case 3: // Bridge
                if (KeyboardState.IsKeyPressed(Keys.Space))
                {
                    // Drop a ball on the bridge
                    var ball = new RigidBody(
                        new ArtemisEngine.Vector2(12 + Random.Shared.NextSingle() * 6, 15),
                        5f, new CircleShape(0.6f));
                    ball.Restitution = 0.4f;
                    World.AddBody(ball);
                }
                break;
        }
    }

    protected override void UpdatePhysics(float deltaTime)
    {
        // Remove bodies that fall too far
        var toRemove = World.Bodies.Where(b => b.Position.Y < -20 && b != _ground).ToList();
        foreach (var body in toRemove)
        {
            World.RemoveBody(body);
        }
    }

    protected override void Render()
    {
        // Draw ground
        if (_ground != null && _ground.Shape is BoxShape groundBox)
        {
            DrawBox(_ground.Position, groundBox.Width, groundBox.Height, _ground.Rotation,
                new Color4(0.3f, 0.5f, 0.3f, 1f));
        }

        // Draw all bodies
        foreach (var body in World.Bodies)
        {
            if (body == _ground) continue;

            Color4 color = GetBodyColor(body);

            if (body.Shape is CircleShape circle)
            {
                DrawCircle(body.Position, circle.Radius, color);
                // Draw rotation indicator
                ArtemisEngine.Vector2 dir = new ArtemisEngine.Vector2(
                    MathF.Cos(body.Rotation),
                    MathF.Sin(body.Rotation)
                );
                DrawLine(body.Position, body.Position + dir * circle.Radius * 0.8f, Color4.White);
            }
            else if (body.Shape is BoxShape box)
            {
                DrawBox(body.Position, box.Width, box.Height, body.Rotation, color);
            }
        }

        // Draw joints
        foreach (var joint in World.Joints)
        {
            DrawJoint(joint);
        }
    }

    private Color4 GetBodyColor(RigidBody body)
    {
        if (body.IsStatic)
            return new Color4(0.4f, 0.4f, 0.5f, 1f);

        // Ragdoll parts
        if (_ragdoll != null)
        {
            if (body == _ragdoll.Head)
                return new Color4(1f, 0.8f, 0.6f, 1f); // Skin color
            if (body == _ragdoll.Torso)
                return new Color4(0.3f, 0.5f, 0.9f, 1f); // Blue shirt
            if (body == _ragdoll.LeftArm || body == _ragdoll.RightArm)
                return new Color4(1f, 0.8f, 0.6f, 1f); // Skin
            if (body == _ragdoll.LeftLeg || body == _ragdoll.RightLeg)
                return new Color4(0.3f, 0.3f, 0.4f, 1f); // Dark pants
        }

        // Vehicle parts
        if (_vehicle != null)
        {
            if (body == _vehicle.Body)
                return new Color4(0.9f, 0.2f, 0.2f, 1f); // Red car
            if (body == _vehicle.FrontWheel || body == _vehicle.BackWheel)
                return new Color4(0.2f, 0.2f, 0.2f, 1f); // Black wheels
        }

        // Chain/bridge
        if (_chain.Contains(body) || _bridge.Contains(body))
            return new Color4(0.6f, 0.4f, 0.2f, 1f); // Brown

        return new Color4(0.7f, 0.7f, 0.7f, 1f); // Default gray
    }

    private void DrawJoint(Joint joint)
    {
        ArtemisEngine.Vector2 posA = joint.BodyA.Position;
        ArtemisEngine.Vector2 posB = joint.BodyB.Position;

        Color4 color = joint switch
        {
            SpringJoint => new Color4(0.2f, 0.9f, 0.2f, 1f),
            DistanceJoint => new Color4(0.2f, 0.8f, 0.9f, 1f),
            RevoluteJoint => new Color4(0.9f, 0.2f, 0.9f, 1f),
            RopeJoint => new Color4(0.9f, 0.9f, 0.2f, 1f),
            _ => new Color4(0.8f, 0.8f, 0.8f, 1f)
        };

        DrawLine(posA, posB, color, 2f);

        // Draw joint anchors as small circles
        DrawCircle(posA, 0.15f, color);
        DrawCircle(posB, 0.15f, color);
    }

    protected override void DrawUI()
    {
        string controls = _currentDemo switch
        {
            0 => "[SPACE]Impulse",
            1 => "[A/D]Drive [SPACE]Jump",
            2 => "[SPACE]Swing",
            3 => "[SPACE]Drop Ball",
            _ => ""
        };

        Title = $"Advanced Physics | Demo: {_demoNames[_currentDemo]} | Bodies: {World.Bodies.Count} | " +
                $"Joints: {World.Joints.Count} | {controls} | " +
                $"[1-4]Demos [R]Reset [P]Pause [ESC]Quit";
    }
}
