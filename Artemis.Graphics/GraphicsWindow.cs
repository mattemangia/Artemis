using OpenTK.Graphics.OpenGL4;
using OpenTK.Mathematics;
using OpenTK.Windowing.Common;
using OpenTK.Windowing.Desktop;
using OpenTK.Windowing.GraphicsLibraryFramework;

namespace Artemis.Graphics;

/// <summary>
/// Base graphics window for physics demos using OpenTK
/// Provides 2D rendering capabilities with camera control
/// </summary>
public abstract class GraphicsWindow : GameWindow
{
    protected PhysicsWorld World { get; set; }
    protected Vector2d CameraPosition { get; set; } = Vector2d.Zero;
    protected double CameraZoom { get; set; } = 20.0;

    private int _shaderProgram;
    private int _vao;
    private int _vbo;

    // Circle rendering
    private int _circleVao;
    private int _circleVbo;
    private const int CircleSegments = 32;

    protected bool IsPaused { get; set; } = false;
    protected float TimeScale { get; set; } = 1.0f;
    protected float FixedTimeStep { get; set; } = 1f / 60f;

    private float _accumulator = 0;

    public GraphicsWindow(int width, int height, string title)
        : base(
            GameWindowSettings.Default,
            new NativeWindowSettings()
            {
                ClientSize = new Vector2i(width, height),
                Title = title,
                APIVersion = new Version(3, 3)
            })
    {
        World = new PhysicsWorld(new Vector2(0, -10f));
    }

    protected override void OnLoad()
    {
        base.OnLoad();

        GL.ClearColor(0.1f, 0.1f, 0.15f, 1.0f);

        // Create simple shader program
        CreateShaders();

        // Create VAO/VBO for line/box rendering
        CreateBuffers();

        // Create circle buffer
        CreateCircleBuffer();

        // Initialize the demo
        Initialize();
    }

    private void CreateShaders()
    {
        string vertexShaderSource = @"
            #version 330 core
            layout (location = 0) in vec2 aPosition;
            uniform mat4 projection;
            uniform vec4 color;
            out vec4 fragColor;
            void main()
            {
                gl_Position = projection * vec4(aPosition, 0.0, 1.0);
                fragColor = color;
            }
        ";

        string fragmentShaderSource = @"
            #version 330 core
            in vec4 fragColor;
            out vec4 FragColor;
            void main()
            {
                FragColor = fragColor;
            }
        ";

        int vertexShader = GL.CreateShader(ShaderType.VertexShader);
        GL.ShaderSource(vertexShader, vertexShaderSource);
        GL.CompileShader(vertexShader);

        int fragmentShader = GL.CreateShader(ShaderType.FragmentShader);
        GL.ShaderSource(fragmentShader, fragmentShaderSource);
        GL.CompileShader(fragmentShader);

        _shaderProgram = GL.CreateProgram();
        GL.AttachShader(_shaderProgram, vertexShader);
        GL.AttachShader(_shaderProgram, fragmentShader);
        GL.LinkProgram(_shaderProgram);

        GL.DeleteShader(vertexShader);
        GL.DeleteShader(fragmentShader);
    }

    private void CreateBuffers()
    {
        _vao = GL.GenVertexArray();
        _vbo = GL.GenBuffer();

        GL.BindVertexArray(_vao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);

        GL.VertexAttribPointer(0, 2, VertexAttribPointerType.Float, false, 2 * sizeof(float), 0);
        GL.EnableVertexAttribArray(0);

        GL.BindVertexArray(0);
    }

    private void CreateCircleBuffer()
    {
        _circleVao = GL.GenVertexArray();
        _circleVbo = GL.GenBuffer();

        // Generate circle vertices
        float[] circleVertices = new float[(CircleSegments + 2) * 2];
        circleVertices[0] = 0; // Center
        circleVertices[1] = 0;

        for (int i = 0; i <= CircleSegments; i++)
        {
            float angle = (float)(2 * Math.PI * i / CircleSegments);
            circleVertices[(i + 1) * 2] = MathF.Cos(angle);
            circleVertices[(i + 1) * 2 + 1] = MathF.Sin(angle);
        }

        GL.BindVertexArray(_circleVao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _circleVbo);
        GL.BufferData(BufferTarget.ArrayBuffer, circleVertices.Length * sizeof(float), circleVertices, BufferUsageHint.StaticDraw);

        GL.VertexAttribPointer(0, 2, VertexAttribPointerType.Float, false, 2 * sizeof(float), 0);
        GL.EnableVertexAttribArray(0);

        GL.BindVertexArray(0);
    }

    protected override void OnUpdateFrame(FrameEventArgs args)
    {
        base.OnUpdateFrame(args);

        if (KeyboardState.IsKeyDown(Keys.Escape))
        {
            Close();
        }

        // Pause toggle
        if (KeyboardState.IsKeyPressed(Keys.P))
        {
            IsPaused = !IsPaused;
        }

        // Time scale controls
        if (KeyboardState.IsKeyDown(Keys.Minus) || KeyboardState.IsKeyDown(Keys.KeyPadSubtract))
        {
            TimeScale = Math.Max(0.1f, TimeScale - 0.01f);
        }
        if (KeyboardState.IsKeyDown(Keys.Equal) || KeyboardState.IsKeyDown(Keys.KeyPadAdd))
        {
            TimeScale = Math.Min(3.0f, TimeScale + 0.01f);
        }

        // Camera controls
        double panSpeed = 0.5;
        if (KeyboardState.IsKeyDown(Keys.Left))
            CameraPosition = new Vector2d(CameraPosition.X - panSpeed, CameraPosition.Y);
        if (KeyboardState.IsKeyDown(Keys.Right))
            CameraPosition = new Vector2d(CameraPosition.X + panSpeed, CameraPosition.Y);
        if (KeyboardState.IsKeyDown(Keys.Up))
            CameraPosition = new Vector2d(CameraPosition.X, CameraPosition.Y + panSpeed);
        if (KeyboardState.IsKeyDown(Keys.Down))
            CameraPosition = new Vector2d(CameraPosition.X, CameraPosition.Y - panSpeed);

        // Zoom controls
        if (KeyboardState.IsKeyDown(Keys.Z))
            CameraZoom *= 1.02;
        if (KeyboardState.IsKeyDown(Keys.X))
            CameraZoom /= 1.02;

        // Handle demo-specific input
        HandleInput();

        // Update physics
        if (!IsPaused)
        {
            _accumulator += (float)args.Time;

            while (_accumulator >= FixedTimeStep)
            {
                World.Step(FixedTimeStep * TimeScale);
                UpdatePhysics(FixedTimeStep * TimeScale);
                _accumulator -= FixedTimeStep;
            }
        }
    }

    protected override void OnRenderFrame(FrameEventArgs args)
    {
        base.OnRenderFrame(args);

        GL.Clear(ClearBufferMask.ColorBufferBit);

        GL.UseProgram(_shaderProgram);

        // Set projection matrix
        Matrix4 projection = CreateProjectionMatrix();
        int projectionLoc = GL.GetUniformLocation(_shaderProgram, "projection");
        GL.UniformMatrix4(projectionLoc, false, ref projection);

        // Render demo content
        Render();

        // Draw UI overlay
        DrawUI();

        SwapBuffers();
    }

    private Matrix4 CreateProjectionMatrix()
    {
        float aspect = (float)Size.X / Size.Y;
        float halfWidth = (float)(CameraZoom * aspect);
        float halfHeight = (float)CameraZoom;

        return Matrix4.CreateOrthographicOffCenter(
            (float)CameraPosition.X - halfWidth,
            (float)CameraPosition.X + halfWidth,
            (float)CameraPosition.Y - halfHeight,
            (float)CameraPosition.Y + halfHeight,
            -1f, 1f);
    }

    protected override void OnResize(ResizeEventArgs e)
    {
        base.OnResize(e);
        GL.Viewport(0, 0, e.Width, e.Height);
    }

    protected override void OnUnload()
    {
        GL.DeleteVertexArray(_vao);
        GL.DeleteBuffer(_vbo);
        GL.DeleteVertexArray(_circleVao);
        GL.DeleteBuffer(_circleVbo);
        GL.DeleteProgram(_shaderProgram);
        base.OnUnload();
    }

    #region Drawing Methods

    protected void DrawCircle(Vector2 center, double radius, Color4 color, bool filled = true)
    {
        int colorLoc = GL.GetUniformLocation(_shaderProgram, "color");
        GL.Uniform4(colorLoc, color);

        // Create transformation for the circle
        float[] vertices = new float[(CircleSegments + 2) * 2];
        float centerX = (float)center.X;
        float centerY = (float)center.Y;
        float radiusF = (float)radius;
        vertices[0] = centerX;
        vertices[1] = centerY;

        for (int i = 0; i <= CircleSegments; i++)
        {
            float angle = (float)(2 * Math.PI * i / CircleSegments);
            vertices[(i + 1) * 2] = centerX + MathF.Cos(angle) * radiusF;
            vertices[(i + 1) * 2 + 1] = centerY + MathF.Sin(angle) * radiusF;
        }

        GL.BindVertexArray(_vao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertices.Length * sizeof(float), vertices, BufferUsageHint.DynamicDraw);

        if (filled)
        {
            GL.DrawArrays(PrimitiveType.TriangleFan, 0, CircleSegments + 2);
        }
        else
        {
            GL.DrawArrays(PrimitiveType.LineLoop, 1, CircleSegments);
        }

        GL.BindVertexArray(0);
    }

    protected void DrawBox(Vector2 center, double width, double height, double rotation, Color4 color, bool filled = true)
    {
        int colorLoc = GL.GetUniformLocation(_shaderProgram, "color");
        GL.Uniform4(colorLoc, color);

        float centerX = (float)center.X;
        float centerY = (float)center.Y;
        float hw = (float)(width / 2);
        float hh = (float)(height / 2);

        // Compute rotated corners
        float cos = MathF.Cos((float)rotation);
        float sin = MathF.Sin((float)rotation);

        Vector2[] corners = new Vector2[4]
        {
            new Vector2(-hw, -hh),
            new Vector2(hw, -hh),
            new Vector2(hw, hh),
            new Vector2(-hw, hh)
        };

        float[] vertices = new float[8];
        for (int i = 0; i < 4; i++)
        {
            float rx = (float)(corners[i].X * cos - corners[i].Y * sin);
            float ry = (float)(corners[i].X * sin + corners[i].Y * cos);
            vertices[i * 2] = centerX + rx;
            vertices[i * 2 + 1] = centerY + ry;
        }

        GL.BindVertexArray(_vao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertices.Length * sizeof(float), vertices, BufferUsageHint.DynamicDraw);

        if (filled)
        {
            GL.DrawArrays(PrimitiveType.TriangleFan, 0, 4);
        }
        else
        {
            GL.DrawArrays(PrimitiveType.LineLoop, 0, 4);
        }

        GL.BindVertexArray(0);
    }

    protected void DrawLine(Vector2 start, Vector2 end, Color4 color, float lineWidth = 1f)
    {
        int colorLoc = GL.GetUniformLocation(_shaderProgram, "color");
        GL.Uniform4(colorLoc, color);

        float[] vertices = new float[]
        {
            (float)start.X, (float)start.Y,
            (float)end.X, (float)end.Y
        };

        GL.LineWidth(lineWidth);

        GL.BindVertexArray(_vao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertices.Length * sizeof(float), vertices, BufferUsageHint.DynamicDraw);

        GL.DrawArrays(PrimitiveType.Lines, 0, 2);

        GL.BindVertexArray(0);
        GL.LineWidth(1f);
    }

    protected void DrawRing(Vector2 center, double radius, Color4 color, float lineWidth = 2f)
    {
        int colorLoc = GL.GetUniformLocation(_shaderProgram, "color");
        GL.Uniform4(colorLoc, color);

        float[] vertices = new float[CircleSegments * 2];
        float centerX = (float)center.X;
        float centerY = (float)center.Y;
        float radiusF = (float)radius;
        for (int i = 0; i < CircleSegments; i++)
        {
            float angle = (float)(2 * Math.PI * i / CircleSegments);
            vertices[i * 2] = centerX + MathF.Cos(angle) * radiusF;
            vertices[i * 2 + 1] = centerY + MathF.Sin(angle) * radiusF;
        }

        GL.LineWidth(lineWidth);

        GL.BindVertexArray(_vao);
        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertices.Length * sizeof(float), vertices, BufferUsageHint.DynamicDraw);

        GL.DrawArrays(PrimitiveType.LineLoop, 0, CircleSegments);

        GL.BindVertexArray(0);
        GL.LineWidth(1f);
    }

    protected void DrawTrail(List<Vector2> points, Color4 startColor, Color4 endColor)
    {
        if (points.Count < 2) return;

        for (int i = 0; i < points.Count - 1; i++)
        {
            float t = (float)i / points.Count;
            Color4 color = new Color4(
                startColor.R * (1 - t) + endColor.R * t,
                startColor.G * (1 - t) + endColor.G * t,
                startColor.B * (1 - t) + endColor.B * t,
                startColor.A * (1 - t) + endColor.A * t
            );

            DrawLine(points[i], points[i + 1], color, 1f);
        }
    }

    protected void DrawBody(RigidBody body, Color4 color)
    {
        if (body.Shape is CircleShape circle)
        {
            DrawCircle(body.Position, circle.Radius, color);

            // Draw rotation indicator
            float rotation = (float)body.Rotation;
            Vector2 dir = new Vector2(
                MathF.Cos(rotation),
                MathF.Sin(rotation)
            );
            DrawLine(body.Position, body.Position + dir * circle.Radius * 0.8f, Color4.White, 1f);
        }
        else if (body.Shape is BoxShape box)
        {
            DrawBox(body.Position, box.Width, box.Height, body.Rotation, color);
        }
    }

    #endregion

    #region Abstract Methods - Implement in derived classes

    protected abstract void Initialize();
    protected abstract void HandleInput();
    protected abstract void UpdatePhysics(float deltaTime);
    protected abstract void Render();
    protected virtual void DrawUI() { }

    #endregion
}
