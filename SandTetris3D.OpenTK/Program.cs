using OpenTK.Graphics.OpenGL4;
using OpenTK.Mathematics;
using OpenTK.Windowing.Common;
using OpenTK.Windowing.Desktop;
using OpenTK.Windowing.GraphicsLibraryFramework;
using Artemis.Demo;
using Artemis.Core;
using Artemis.Physics2D;

namespace SandTetris3D.OpenTK;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== Sand Tetris 3D - OpenTK Edition (Block Mode) ===");
        Console.WriteLine("A 3D physics-based puzzle game with Blocks");

        var settings = new NativeWindowSettings()
        {
            ClientSize = new Vector2i(1280, 720),
            Title = "Sand Tetris 3D - Artemis Physics Engine",
            APIVersion = new Version(3, 3)
        };

        using var window = new SandTetrisWindow(GameWindowSettings.Default, settings);
        window.Run();
    }
}

public class SandTetrisWindow : GameWindow
{
    private Artemis.Demo.SandTetris3D _game = null!;

    // Camera
    private Vector3 _cameraPosition = new Vector3(12f, 15f, 18f);
    private Vector3 _cameraTarget = new Vector3(0, 7.5f, 0);
    private float _cameraRotation = 0f;

    // Drop position
    private double _dropX = 0;
    private double _dropZ = 0;
    private bool _waitingForSettle = false;
    private int _settleFrames = 0;
    private const int SettleFramesRequired = 60;

    // Mouse control
    private Vector2 _lastMousePos;
    private bool _isRotating = false;

    // Shaders
    private int _shaderProgram;
    private int _cubeVao;
    private int _cubeVbo;

    // Physics accumulator
    private double _physicsAccumulator = 0;
    private const double PhysicsStep = 1.0 / 60.0;

    public SandTetrisWindow(GameWindowSettings gameWindowSettings, NativeWindowSettings nativeWindowSettings)
        : base(gameWindowSettings, nativeWindowSettings)
    {
    }

    protected override void OnLoad()
    {
        base.OnLoad();

        GL.ClearColor(0.05f, 0.08f, 0.12f, 1.0f);
        GL.Enable(EnableCap.DepthTest);
        GL.Enable(EnableCap.CullFace);

        // Initialize game
        _game = new Artemis.Demo.SandTetris3D(true);

        // Create shaders
        CreateShaders();

        // Create cube geometry (reused for walls and blocks)
        CreateCube();
    }

    private void CreateShaders()
    {
        string vertexSource = @"
            #version 330 core
            layout(location = 0) in vec3 aPosition;
            layout(location = 1) in vec3 aNormal;

            uniform mat4 model;
            uniform mat4 view;
            uniform mat4 projection;

            out vec3 FragPos;
            out vec3 Normal;

            void main()
            {
                FragPos = vec3(model * vec4(aPosition, 1.0));
                Normal = mat3(transpose(inverse(model))) * aNormal;
                gl_Position = projection * view * model * vec4(aPosition, 1.0);
            }
        ";

        string fragmentSource = @"
            #version 330 core
            in vec3 FragPos;
            in vec3 Normal;

            uniform vec3 lightPos;
            uniform vec3 viewPos;
            uniform vec3 objectColor;
            uniform float shineProgress;

            out vec4 FragColor;

            void main()
            {
                // Ambient
                float ambientStrength = 0.3;
                vec3 ambient = ambientStrength * vec3(1.0);

                // Diffuse
                vec3 norm = normalize(Normal);
                vec3 lightDir = normalize(lightPos - FragPos);
                float diff = max(dot(norm, lightDir), 0.0);
                vec3 diffuse = diff * vec3(1.0);

                // Specular
                float specularStrength = 0.5;
                vec3 viewDir = normalize(viewPos - FragPos);
                vec3 reflectDir = reflect(-lightDir, norm);
                float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
                vec3 specular = specularStrength * spec * vec3(1.0);

                // Shine effect for particles being removed
                vec3 color = objectColor;
                if (shineProgress > 0.0) {
                    float pulse = sin(shineProgress * 3.14159 * 4.0) * 0.5 + 0.5;
                    color = mix(objectColor, vec3(1.0), pulse * 0.5);
                }

                vec3 result = (ambient + diffuse + specular) * color;
                FragColor = vec4(result, 1.0);
            }
        ";

        int vertexShader = GL.CreateShader(ShaderType.VertexShader);
        GL.ShaderSource(vertexShader, vertexSource);
        GL.CompileShader(vertexShader);

        int fragmentShader = GL.CreateShader(ShaderType.FragmentShader);
        GL.ShaderSource(fragmentShader, fragmentSource);
        GL.CompileShader(fragmentShader);

        _shaderProgram = GL.CreateProgram();
        GL.AttachShader(_shaderProgram, vertexShader);
        GL.AttachShader(_shaderProgram, fragmentShader);
        GL.LinkProgram(_shaderProgram);

        GL.DeleteShader(vertexShader);
        GL.DeleteShader(fragmentShader);
    }

    private void CreateCube()
    {
        float[] vertices = {
            // Positions          // Normals
            // Front face (CCW)
            -0.5f, -0.5f,  0.5f,  0f, 0f, 1f,
             0.5f, -0.5f,  0.5f,  0f, 0f, 1f,
             0.5f,  0.5f,  0.5f,  0f, 0f, 1f,
            -0.5f,  0.5f,  0.5f,  0f, 0f, 1f,
            // Back face (Corrected to CCW)
            -0.5f, -0.5f, -0.5f,  0f, 0f, -1f,
            -0.5f,  0.5f, -0.5f,  0f, 0f, -1f,
             0.5f,  0.5f, -0.5f,  0f, 0f, -1f,
             0.5f, -0.5f, -0.5f,  0f, 0f, -1f,
            // Left face (CCW)
            -0.5f,  0.5f,  0.5f, -1f, 0f, 0f,
            -0.5f,  0.5f, -0.5f, -1f, 0f, 0f,
            -0.5f, -0.5f, -0.5f, -1f, 0f, 0f,
            -0.5f, -0.5f,  0.5f, -1f, 0f, 0f,
            // Right face (Corrected to CCW)
             0.5f,  0.5f,  0.5f,  1f, 0f, 0f,
             0.5f, -0.5f,  0.5f,  1f, 0f, 0f,
             0.5f, -0.5f, -0.5f,  1f, 0f, 0f,
             0.5f,  0.5f, -0.5f,  1f, 0f, 0f,
             // Top face (Corrected to CCW)
            -0.5f,  0.5f, -0.5f,  0f, 1f, 0f,
            -0.5f,  0.5f,  0.5f,  0f, 1f, 0f,
             0.5f,  0.5f,  0.5f,  0f, 1f, 0f,
             0.5f,  0.5f, -0.5f,  0f, 1f, 0f,
             // Bottom face (CCW)
            -0.5f, -0.5f, -0.5f,  0f, -1f, 0f,
             0.5f, -0.5f, -0.5f,  0f, -1f, 0f,
             0.5f, -0.5f,  0.5f,  0f, -1f, 0f,
            -0.5f, -0.5f,  0.5f,  0f, -1f, 0f,
        };

        _cubeVao = GL.GenVertexArray();
        _cubeVbo = GL.GenBuffer();

        GL.BindVertexArray(_cubeVao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _cubeVbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertices.Length * sizeof(float), vertices, BufferUsageHint.StaticDraw);

        GL.VertexAttribPointer(0, 3, VertexAttribPointerType.Float, false, 6 * sizeof(float), 0);
        GL.EnableVertexAttribArray(0);

        GL.VertexAttribPointer(1, 3, VertexAttribPointerType.Float, false, 6 * sizeof(float), 3 * sizeof(float));
        GL.EnableVertexAttribArray(1);

        GL.BindVertexArray(0);
    }

    protected override void OnUpdateFrame(FrameEventArgs args)
    {
        base.OnUpdateFrame(args);

        if (KeyboardState.IsKeyDown(Keys.Escape))
            Close();

        // Mouse-based camera rotation (right mouse button)
        var mousePos = MouseState.Position;
        if (MouseState.IsButtonDown(MouseButton.Right))
        {
            if (!_isRotating)
            {
                _isRotating = true;
                _lastMousePos = mousePos;
            }
            else
            {
                float deltaX = mousePos.X - _lastMousePos.X;
                _cameraRotation -= deltaX * 0.01f;
                _lastMousePos = mousePos;
            }
        }
        else
        {
            _isRotating = false;
        }

        // Keyboard camera rotation (fallback)
        if (KeyboardState.IsKeyDown(Keys.Q))
            _cameraRotation += 0.02f;
        if (KeyboardState.IsKeyDown(Keys.E))
            _cameraRotation -= 0.02f;

        // Update camera position based on rotation
        float radius = 22f;
        _cameraPosition = new Vector3(
            MathF.Sin(_cameraRotation) * radius,
            15f,
            MathF.Cos(_cameraRotation) * radius
        );

        // Mouse-based aiming (move mouse to aim, left click to drop)
        if (!_waitingForSettle && !_game.IsGameOver)
        {
            // Convert mouse position to drop coordinates
            float normalizedX = mousePos.X / Size.X;  // 0 to 1
            float normalizedY = mousePos.Y / Size.Y;  // 0 to 1

            _dropX = (normalizedX - 0.5) * 9.0;  // -4.5 to 4.5
            _dropZ = 0; // Fixed Z for 2D logic
            _dropX = Math.Clamp(_dropX, -4.5, 4.5);

            // Left click to drop ball
            if (MouseState.IsButtonPressed(MouseButton.Left))
            {
                _game.DropBall(_dropX, _dropZ);
                _waitingForSettle = true;
                _settleFrames = 0;
            }
        }

        // Keyboard aiming controls (fallback)
        if (!_waitingForSettle)
        {
            double moveStep = 0.3;
            if (KeyboardState.IsKeyDown(Keys.Left))
                _dropX = Math.Max(_dropX - moveStep, -4.5);
            if (KeyboardState.IsKeyDown(Keys.Right))
                _dropX = Math.Min(_dropX + moveStep, 4.5);

            // Space to drop ball (keyboard fallback)
            if (KeyboardState.IsKeyPressed(Keys.Space) && !_game.IsGameOver)
            {
                _game.DropBall(_dropX, _dropZ);
                _waitingForSettle = true;
                _settleFrames = 0;
            }
        }

        // Reset
        if (KeyboardState.IsKeyPressed(Keys.R))
        {
            _game.Reset();
            _dropX = 0;
            _dropZ = 0;
            _waitingForSettle = false;
        }

        // Physics update - ALWAYS run physics simulation (not just when settling)
        _physicsAccumulator += args.Time;
        while (_physicsAccumulator >= PhysicsStep)
        {
            _game.Update(PhysicsStep);
            _physicsAccumulator -= PhysicsStep;

            if (_waitingForSettle)
                _settleFrames++;
        }

        // Check if settled
        if (_waitingForSettle && _settleFrames >= SettleFramesRequired)
        {
            _game.CheckAndClearLines();
            _game.CheckGameState();

            if (!_game.IsGameOver)
            {
                _game.NextTurn();
            }

            _waitingForSettle = false;
            _settleFrames = 0;
        }

        // Update title
        UpdateTitle();
    }

    private void UpdateTitle()
    {
        string status = _game.IsGameOver
            ? (_game.IsVictory ? " [VICTORY!]" : " [GAME OVER]")
            : (_waitingForSettle ? " [Settling...]" : "");

        Title = $"Sand Tetris 3D (Blocks) | Turn: {_game.Turn} | Score: {_game.Score} | " +
                $"Blocks: {_game.ActiveParticleCount} | " +
                $"Ball: {_game.CurrentBall.ColorName} ({_game.CurrentBall.Radius:F1}){status} | " +
                "[Mouse]Aim [LClick]Drop [RDrag]Rotate [R]Reset";
    }

    protected override void OnRenderFrame(FrameEventArgs args)
    {
        base.OnRenderFrame(args);

        GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);

        GL.UseProgram(_shaderProgram);

        // Set up matrices
        Matrix4 view = Matrix4.LookAt(_cameraPosition, _cameraTarget, Vector3.UnitY);
        Matrix4 projection = Matrix4.CreatePerspectiveFieldOfView(
            MathHelper.DegreesToRadians(45f),
            (float)Size.X / Size.Y,
            0.1f, 100f);

        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "view"), false, ref view);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "projection"), false, ref projection);

        // Light and camera position
        Vector3 lightPos = new Vector3(10, 20, 10);
        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "lightPos"), lightPos);
        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "viewPos"), _cameraPosition);

        // Draw terrarium
        DrawTerrarium();

        // Draw particles/blocks
        DrawBlocks();

        // Draw drop indicator
        if (!_waitingForSettle && !_game.IsGameOver)
        {
            DrawDropIndicator();
        }

        SwapBuffers();
    }

    private void DrawTerrarium()
    {
        var (width, height, depth) = _game.GetTerrariumSize();

        // Draw floor
        Matrix4 floorModel = Matrix4.CreateScale((float)width, 0.1f, (float)depth);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "model"), false, ref floorModel);
        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "objectColor"), 0.2f, 0.3f, 0.2f);
        GL.Uniform1(GL.GetUniformLocation(_shaderProgram, "shineProgress"), 0f);

        GL.BindVertexArray(_cubeVao);
        GL.DrawArrays(PrimitiveType.TriangleFan, 0, 4); // Top
        GL.DrawArrays(PrimitiveType.TriangleFan, 16, 4); // Bottom? Cube has 24 vertices

        // Draw walls (wireframe style)
        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line);
        GL.LineWidth(2f);

        Matrix4 wallsModel = Matrix4.CreateScale((float)width, (float)height, (float)depth) *
                             Matrix4.CreateTranslation(0, (float)height / 2, 0);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "model"), false, ref wallsModel);
        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "objectColor"), 0.4f, 0.5f, 0.6f);

        // Draw 12 lines of the cube
        GL.DrawArrays(PrimitiveType.Lines, 0, 24); // Not efficient but works if lines are defined
        // Actually CreateCube uses TRIANGLES (via DrawArrays or ElementBuffer?)
        // The CreateCube above uses raw vertices without indices. 6 faces * 4 vertices = 24.
        // Wait, CreateCube above has 24 vertices (4 per face, 6 faces).
        // GL_TRIANGLE_FAN for each face? No, CreateCube didn't set up Triangle Fans correctly for all faces in one go.
        // It's 4 vertices per face.

        // Re-draw as lines:
        for (int i = 0; i < 6; i++)
        {
            GL.DrawArrays(PrimitiveType.LineLoop, i * 4, 4);
        }

        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
    }

    private void DrawBlocks()
    {
        GL.BindVertexArray(_cubeVao);

        int shineProgressLoc = GL.GetUniformLocation(_shaderProgram, "shineProgress");
        int colorLoc = GL.GetUniformLocation(_shaderProgram, "objectColor");
        int modelLoc = GL.GetUniformLocation(_shaderProgram, "model");

        var bodies = _game.World.Bodies;
        float depth = (float)_game.GetTerrariumSize().depth;

        for (int i = 0; i < bodies.Count; i++)
        {
            var body = bodies[i];
            if (!body.IsActive || body.BodyType == BodyType2D.Static) continue;

            uint color = (body.UserData is uint c) ? c : 0xFFFFFFFF;

            // Check shine
            float shineProgress = 0f;
            if (_game.ShiningBodies.TryGetValue(body.Id, out float progress))
            {
                shineProgress = progress;
                color = _game.GetShineColor(color, progress);
            }

            var (r, g, b) = Artemis.Demo.SandTetris3D.GetRGB(color);

            GL.Uniform1(shineProgressLoc, shineProgress);
            GL.Uniform3(colorLoc, r / 255f, g / 255f, b / 255f);

            // Get dimensions from shape
            float width = 1.0f, height = 1.0f;
            if (body.Shape is Artemis.Physics2D.BoxShape box)
            {
                width = (float)box.Width;
                height = (float)box.Height;
            }

            Matrix4 model = Matrix4.CreateScale(width, height, depth) *
                            Matrix4.CreateRotationZ((float)body.Rotation) *
                            Matrix4.CreateTranslation(
                               (float)body.Position.X,
                               (float)body.Position.Y,
                               0); // Z is 0 centered

            GL.UniformMatrix4(modelLoc, false, ref model);

            // Draw cube (6 faces, 4 vertices each = 24 vertices)
            // Using TriangleFan for each face
            for (int f = 0; f < 6; f++)
            {
                GL.DrawArrays(PrimitiveType.TriangleFan, f * 4, 4);
            }
        }
    }

    private void DrawDropIndicator()
    {
        // Draw a wireframe sphere/box at drop position
        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line);

        var ball = _game.CurrentBall;
        var (r, g, b) = Artemis.Demo.SandTetris3D.GetRGB(ball.Color);

        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "objectColor"), r / 255f, g / 255f, b / 255f);
        GL.Uniform1(GL.GetUniformLocation(_shaderProgram, "shineProgress"), 0f);

        // Represent the drop area
        float radius = (float)ball.Radius;
        float depth = (float)_game.GetTerrariumSize().depth;

        Matrix4 model = Matrix4.CreateScale(radius * 2, radius * 2, depth) *
                       Matrix4.CreateTranslation((float)_dropX, 14f, 0);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "model"), false, ref model);

        GL.BindVertexArray(_cubeVao);
        for (int f = 0; f < 6; f++)
        {
            GL.DrawArrays(PrimitiveType.LineLoop, f * 4, 4);
        }

        // Draw drop line
        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
    }

    protected override void OnResize(ResizeEventArgs e)
    {
        base.OnResize(e);
        GL.Viewport(0, 0, e.Width, e.Height);
    }

    protected override void OnUnload()
    {
        GL.DeleteVertexArray(_cubeVao);
        GL.DeleteBuffer(_cubeVbo);
        GL.DeleteProgram(_shaderProgram);
        base.OnUnload();
    }
}
