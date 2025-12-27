using OpenTK.Graphics.OpenGL4;
using OpenTK.Mathematics;
using OpenTK.Windowing.Common;
using OpenTK.Windowing.Desktop;
using OpenTK.Windowing.GraphicsLibraryFramework;
using Artemis.Demo;
using Artemis.Core;

namespace SandTetris3D.OpenTK;

class Program
{
    static void Main(string[] args)
    {
        Console.WriteLine("=== Sand Tetris 3D - OpenTK Edition ===");
        Console.WriteLine("A 3D physics-based puzzle game");

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
    private int _sphereVao;
    private int _sphereVbo;
    private int _sphereEbo;
    private int _sphereIndexCount;

    // Cube for terrarium walls
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

        // Create sphere geometry for particles
        CreateSphere();

        // Create cube for terrarium walls
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

    private void CreateSphere()
    {
        const int stacks = 12;
        const int slices = 16;

        List<float> vertices = new();
        List<uint> indices = new();

        for (int i = 0; i <= stacks; i++)
        {
            float phi = MathF.PI * i / stacks;
            for (int j = 0; j <= slices; j++)
            {
                float theta = 2 * MathF.PI * j / slices;

                float x = MathF.Sin(phi) * MathF.Cos(theta);
                float y = MathF.Cos(phi);
                float z = MathF.Sin(phi) * MathF.Sin(theta);

                // Position
                vertices.Add(x);
                vertices.Add(y);
                vertices.Add(z);

                // Normal (same as position for unit sphere)
                vertices.Add(x);
                vertices.Add(y);
                vertices.Add(z);
            }
        }

        for (int i = 0; i < stacks; i++)
        {
            for (int j = 0; j < slices; j++)
            {
                uint a = (uint)(i * (slices + 1) + j);
                uint b = (uint)(a + slices + 1);

                indices.Add(a);
                indices.Add(b);
                indices.Add(a + 1);

                indices.Add(b);
                indices.Add(b + 1);
                indices.Add(a + 1);
            }
        }

        _sphereVao = GL.GenVertexArray();
        _sphereVbo = GL.GenBuffer();
        _sphereEbo = GL.GenBuffer();

        GL.BindVertexArray(_sphereVao);

        GL.BindBuffer(BufferTarget.ArrayBuffer, _sphereVbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertices.Count * sizeof(float),
            vertices.ToArray(), BufferUsageHint.StaticDraw);

        GL.BindBuffer(BufferTarget.ElementArrayBuffer, _sphereEbo);
        GL.BufferData(BufferTarget.ElementArrayBuffer, indices.Count * sizeof(uint),
            indices.ToArray(), BufferUsageHint.StaticDraw);

        GL.VertexAttribPointer(0, 3, VertexAttribPointerType.Float, false, 6 * sizeof(float), 0);
        GL.EnableVertexAttribArray(0);

        GL.VertexAttribPointer(1, 3, VertexAttribPointerType.Float, false, 6 * sizeof(float), 3 * sizeof(float));
        GL.EnableVertexAttribArray(1);

        _sphereIndexCount = indices.Count;

        GL.BindVertexArray(0);
    }

    private void CreateCube()
    {
        float[] vertices = {
            // Positions          // Normals
            // Front face
            -0.5f, -0.5f,  0.5f,  0f, 0f, 1f,
             0.5f, -0.5f,  0.5f,  0f, 0f, 1f,
             0.5f,  0.5f,  0.5f,  0f, 0f, 1f,
            -0.5f,  0.5f,  0.5f,  0f, 0f, 1f,
            // Back face
            -0.5f, -0.5f, -0.5f,  0f, 0f, -1f,
             0.5f, -0.5f, -0.5f,  0f, 0f, -1f,
             0.5f,  0.5f, -0.5f,  0f, 0f, -1f,
            -0.5f,  0.5f, -0.5f,  0f, 0f, -1f,
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
            // Map screen X (0 to Width) to drop X (-4.5 to 4.5)
            // Map screen Y (0 to Height) to drop Z (-0.8 to 0.8)
            float normalizedX = mousePos.X / Size.X;  // 0 to 1
            float normalizedY = mousePos.Y / Size.Y;  // 0 to 1

            _dropX = (normalizedX - 0.5) * 9.0;  // -4.5 to 4.5
            _dropZ = (normalizedY - 0.5) * 1.6;  // -0.8 to 0.8
            _dropX = Math.Clamp(_dropX, -4.5, 4.5);
            _dropZ = Math.Clamp(_dropZ, -0.8, 0.8);

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
            if (KeyboardState.IsKeyDown(Keys.Up))
                _dropZ = Math.Max(_dropZ - moveStep, -0.8);
            if (KeyboardState.IsKeyDown(Keys.Down))
                _dropZ = Math.Min(_dropZ + moveStep, 0.8);

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

        Title = $"Sand Tetris 3D | Turn: {_game.Turn} | Score: {_game.Score} | " +
                $"Grains: {_game.ActiveParticleCount} | " +
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

        // Draw particles
        DrawParticles();

        // Draw drop indicator
        if (!_waitingForSettle && !_game.IsGameOver)
        {
            DrawDropIndicator();
        }

        // Draw next ball preview
        DrawNextBallPreview();

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
        GL.DrawArrays(PrimitiveType.TriangleFan, 0, 4);
        GL.DrawArrays(PrimitiveType.TriangleFan, 4, 4);

        // Draw walls (wireframe style - just outline)
        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line);
        GL.LineWidth(2f);

        Matrix4 wallsModel = Matrix4.CreateScale((float)width, (float)height, (float)depth) *
                             Matrix4.CreateTranslation(0, (float)height / 2, 0);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "model"), false, ref wallsModel);
        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "objectColor"), 0.4f, 0.5f, 0.6f);

        GL.BindVertexArray(_cubeVao);
        GL.DrawArrays(PrimitiveType.LineLoop, 0, 4);
        GL.DrawArrays(PrimitiveType.LineLoop, 4, 4);

        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
    }

    private void DrawParticles()
    {
        GL.BindVertexArray(_sphereVao);

        int shineProgressLoc = GL.GetUniformLocation(_shaderProgram, "shineProgress");
        int colorLoc = GL.GetUniformLocation(_shaderProgram, "objectColor");
        int modelLoc = GL.GetUniformLocation(_shaderProgram, "model");

        var particles = _game.Simulation.Particles;
        float radius = (float)_game.Simulation.ParticleRadius;

        for (int i = 0; i < particles.Count; i++)
        {
            var particle = particles[i];
            if (!particle.IsActive) continue;

            // Get color
            var (r, g, b) = Artemis.Demo.SandTetris3D.GetRGB(particle.Color);
            uint color = particle.Color;

            // Check if shining
            float shineProgress = 0f;
            if (_game.ShiningParticles.TryGetValue(i, out float progress))
            {
                shineProgress = progress;
                color = _game.GetShineColor(i, particle.Color);
                (r, g, b) = Artemis.Demo.SandTetris3D.GetRGB(color);
            }

            GL.Uniform1(shineProgressLoc, shineProgress);
            GL.Uniform3(colorLoc, r / 255f, g / 255f, b / 255f);

            Matrix4 model = Matrix4.CreateScale(radius) *
                           Matrix4.CreateTranslation(
                               (float)particle.Position.X,
                               (float)particle.Position.Y,
                               (float)particle.Position.Z);
            GL.UniformMatrix4(modelLoc, false, ref model);

            GL.DrawElements(PrimitiveType.Triangles, _sphereIndexCount, DrawElementsType.UnsignedInt, 0);
        }
    }

    private void DrawDropIndicator()
    {
        // Draw a wireframe sphere at drop position
        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line);

        var ball = _game.CurrentBall;
        var (r, g, b) = Artemis.Demo.SandTetris3D.GetRGB(ball.Color);

        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "objectColor"), r / 255f, g / 255f, b / 255f);
        GL.Uniform1(GL.GetUniformLocation(_shaderProgram, "shineProgress"), 0f);

        Matrix4 model = Matrix4.CreateScale((float)ball.Radius) *
                       Matrix4.CreateTranslation((float)_dropX, 14f, (float)_dropZ);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "model"), false, ref model);

        GL.BindVertexArray(_sphereVao);
        GL.DrawElements(PrimitiveType.Triangles, _sphereIndexCount, DrawElementsType.UnsignedInt, 0);

        // Draw drop line
        GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
    }

    private void DrawNextBallPreview()
    {
        // Draw the next ball floating above the terrarium
        var ball = _game.CurrentBall;
        var (r, g, b) = Artemis.Demo.SandTetris3D.GetRGB(ball.Color);

        GL.Uniform3(GL.GetUniformLocation(_shaderProgram, "objectColor"), r / 255f, g / 255f, b / 255f);
        GL.Uniform1(GL.GetUniformLocation(_shaderProgram, "shineProgress"), 0f);

        // Animate the preview ball
        float bob = MathF.Sin((float)GLFW.GetTime() * 2f) * 0.3f;

        Matrix4 model = Matrix4.CreateScale((float)ball.Radius * 0.5f) *
                       Matrix4.CreateTranslation(-6f, 12f + bob, 0);
        GL.UniformMatrix4(GL.GetUniformLocation(_shaderProgram, "model"), false, ref model);

        GL.BindVertexArray(_sphereVao);
        GL.DrawElements(PrimitiveType.Triangles, _sphereIndexCount, DrawElementsType.UnsignedInt, 0);
    }

    protected override void OnResize(ResizeEventArgs e)
    {
        base.OnResize(e);
        GL.Viewport(0, 0, e.Width, e.Height);
    }

    protected override void OnUnload()
    {
        GL.DeleteVertexArray(_sphereVao);
        GL.DeleteBuffer(_sphereVbo);
        GL.DeleteBuffer(_sphereEbo);
        GL.DeleteVertexArray(_cubeVao);
        GL.DeleteBuffer(_cubeVbo);
        GL.DeleteProgram(_shaderProgram);
        base.OnUnload();
    }
}
