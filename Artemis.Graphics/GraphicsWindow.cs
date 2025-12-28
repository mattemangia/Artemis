using Artemis.Physics2D;
using FontStashSharp;
using OpenTK.Graphics.OpenGL4;
using OpenTK.Mathematics;
using OpenTK.Windowing.Common;
using OpenTK.Windowing.Desktop;
using OpenTK.Windowing.GraphicsLibraryFramework;

namespace Artemis.Graphics;

/// <summary>
/// Base graphics window for physics demos using OpenTK
/// Provides 2D rendering capabilities with camera control and text rendering
/// </summary>
public abstract class GraphicsWindow : GameWindow
{
    protected PhysicsWorld World { get; set; }
    protected Vector2d CameraPosition { get; set; } = Vector2d.Zero;
    protected double CameraZoom { get; set; } = 20.0;

    protected int ShaderProgram { get; private set; }
    private int _vao;
    private int _vbo;

    // Circle rendering
    private int _circleVao;
    private int _circleVbo;
    private const int CircleSegments = 32;

    // Text rendering
    private TextRenderer? _textRenderer;
    protected bool TextRenderingAvailable => _textRenderer?.IsAvailable ?? false;

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

        // Initialize text renderer
        try
        {
            _textRenderer = new TextRenderer(Size.X, Size.Y, 16);
            if (_textRenderer.IsAvailable)
            {
                Console.WriteLine("Text rendering initialized successfully");
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Text rendering unavailable: {ex.Message}");
            _textRenderer = null;
        }

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

        ShaderProgram = GL.CreateProgram();
        GL.AttachShader(ShaderProgram, vertexShader);
        GL.AttachShader(ShaderProgram, fragmentShader);
        GL.LinkProgram(ShaderProgram);

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

        GL.UseProgram(ShaderProgram);

        // Set projection matrix
        Matrix4 projection = CreateProjectionMatrix();
        int projectionLoc = GL.GetUniformLocation(ShaderProgram, "projection");
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
        _textRenderer?.UpdateWindowSize(e.Width, e.Height);
    }

    protected override void OnUnload()
    {
        _textRenderer?.Dispose();
        GL.DeleteVertexArray(_vao);
        GL.DeleteBuffer(_vbo);
        GL.DeleteVertexArray(_circleVao);
        GL.DeleteBuffer(_circleVbo);
        GL.DeleteProgram(ShaderProgram);
        base.OnUnload();
    }

    #region Drawing Methods

    protected void DrawCircle(Vector2 center, double radius, Color4 color, bool filled = true)
    {
        int colorLoc = GL.GetUniformLocation(ShaderProgram, "color");
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
        int colorLoc = GL.GetUniformLocation(ShaderProgram, "color");
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
        int colorLoc = GL.GetUniformLocation(ShaderProgram, "color");
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

    // Helper to convert Vector2D to Vector2
    protected static Vector2 ToVec2(Vector2D v) => new Vector2((float)v.X, (float)v.Y);

    protected void DrawRing(Vector2 center, double radius, Color4 color, float lineWidth = 2f)
    {
        int colorLoc = GL.GetUniformLocation(ShaderProgram, "color");
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

    #region Text Drawing Methods

    /// <summary>
    /// Draw text at screen coordinates (pixels from top-left)
    /// </summary>
    protected void DrawScreenText(string text, float x, float y, Color4 color)
    {
        if (_textRenderer == null || !_textRenderer.IsAvailable) return;
        _textRenderer.DrawText(text, x, y, ToFSColor(color));
    }

    /// <summary>
    /// Draw text at screen coordinates with custom font size
    /// </summary>
    protected void DrawScreenText(string text, float x, float y, int fontSize, Color4 color)
    {
        if (_textRenderer == null || !_textRenderer.IsAvailable) return;
        _textRenderer.DrawText(text, x, y, fontSize, ToFSColor(color));
    }

    /// <summary>
    /// Draw text at world coordinates (will be transformed by camera)
    /// </summary>
    protected void DrawWorldText(string text, Vector2 worldPos, Color4 color)
    {
        if (_textRenderer == null || !_textRenderer.IsAvailable) return;
        var screenPos = WorldToScreen(worldPos);
        _textRenderer.DrawText(text, screenPos.X, screenPos.Y, ToFSColor(color));
    }

    /// <summary>
    /// Draw text at world coordinates with custom font size
    /// </summary>
    protected void DrawWorldText(string text, Vector2 worldPos, int fontSize, Color4 color)
    {
        if (_textRenderer == null || !_textRenderer.IsAvailable) return;
        var screenPos = WorldToScreen(worldPos);
        _textRenderer.DrawText(text, screenPos.X, screenPos.Y, fontSize, ToFSColor(color));
    }

    /// <summary>
    /// Draw centered text at screen coordinates
    /// </summary>
    protected void DrawScreenTextCentered(string text, float x, float y, int fontSize, Color4 color)
    {
        if (_textRenderer == null || !_textRenderer.IsAvailable) return;
        var size = _textRenderer.MeasureText(text, fontSize);
        _textRenderer.DrawText(text, x - size.X / 2, y - size.Y / 2, fontSize, ToFSColor(color));
    }

    /// <summary>
    /// Convert world coordinates to screen coordinates
    /// </summary>
    protected System.Numerics.Vector2 WorldToScreen(Vector2 worldPos)
    {
        float aspect = (float)Size.X / Size.Y;
        float halfWidth = (float)(CameraZoom * aspect);
        float halfHeight = (float)CameraZoom;

        // Normalized coordinates (-1 to 1)
        float nx = (float)((worldPos.X - CameraPosition.X) / halfWidth);
        float ny = (float)((worldPos.Y - CameraPosition.Y) / halfHeight);

        // Screen coordinates
        float sx = (nx + 1) * 0.5f * Size.X;
        float sy = (1 - ny) * 0.5f * Size.Y; // Flip Y for screen space

        return new System.Numerics.Vector2(sx, sy);
    }

    /// <summary>
    /// Measure text size in pixels
    /// </summary>
    protected System.Numerics.Vector2 MeasureText(string text, int fontSize = 16)
    {
        if (_textRenderer == null || !_textRenderer.IsAvailable)
            return System.Numerics.Vector2.Zero;
        return _textRenderer.MeasureText(text, fontSize);
    }

    private static FSColor ToFSColor(Color4 color)
    {
        return new FSColor(
            (byte)(color.R * 255),
            (byte)(color.G * 255),
            (byte)(color.B * 255),
            (byte)(color.A * 255));
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
