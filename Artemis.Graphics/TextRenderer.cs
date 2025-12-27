using FontStashSharp;
using FontStashSharp.Interfaces;
using OpenTK.Graphics.OpenGL4;
using OpenTK.Mathematics;
using System.Numerics;

namespace Artemis.Graphics;

/// <summary>
/// OpenGL text renderer using FontStashSharp
/// </summary>
public class TextRenderer : IDisposable
{
    private readonly FontSystem _fontSystem;
    private readonly SpriteFontBase _font;
    private readonly TextureRenderer _textureRenderer;
    private int _windowWidth;
    private int _windowHeight;

    public TextRenderer(int windowWidth, int windowHeight, int fontSize = 16)
    {
        _windowWidth = windowWidth;
        _windowHeight = windowHeight;

        var settings = new FontSystemSettings
        {
            FontResolutionFactor = 2,
            KernelWidth = 2,
            KernelHeight = 2
        };

        _fontSystem = new FontSystem(settings);

        // Use built-in default font data
        byte[] fontData = GetEmbeddedFontData();
        _fontSystem.AddFont(fontData);

        _font = _fontSystem.GetFont(fontSize);
        _textureRenderer = new TextureRenderer();
    }

    private byte[] GetEmbeddedFontData()
    {
        // Generate a simple bitmap font (basic ASCII characters)
        // For production, you would embed a TTF file
        // Using a minimal embedded font approach
        return GenerateDefaultFontData();
    }

    private byte[] GenerateDefaultFontData()
    {
        // This creates a minimal valid TTF font structure
        // In practice, you'd bundle a proper TTF file like Roboto or DejaVu Sans
        // For now, we'll use the system fallback

        // Try to load from common system font paths
        string[] fontPaths = new[]
        {
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/TTF/DejaVuSans.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
            "C:\\Windows\\Fonts\\arial.ttf",
            "C:\\Windows\\Fonts\\segoeui.ttf",
            "/System/Library/Fonts/Helvetica.ttc"
        };

        foreach (var path in fontPaths)
        {
            if (File.Exists(path))
            {
                return File.ReadAllBytes(path);
            }
        }

        // Fallback: return minimal font data (will use built-in fallback)
        throw new FileNotFoundException("No suitable font found. Please install DejaVu or Liberation fonts.");
    }

    public void UpdateWindowSize(int width, int height)
    {
        _windowWidth = width;
        _windowHeight = height;
    }

    public void DrawText(string text, float x, float y, FSColor color)
    {
        _textureRenderer.Begin(_windowWidth, _windowHeight);
        _font.DrawText(_textureRenderer, text, new System.Numerics.Vector2(x, y), color);
        _textureRenderer.End();
    }

    public void DrawText(string text, float x, float y, float scale, FSColor color)
    {
        var scaledFont = _fontSystem.GetFont((int)(_font.FontSize * scale));
        _textureRenderer.Begin(_windowWidth, _windowHeight);
        scaledFont.DrawText(_textureRenderer, text, new System.Numerics.Vector2(x, y), color);
        _textureRenderer.End();
    }

    public void Dispose()
    {
        _fontSystem.Dispose();
        _textureRenderer.Dispose();
    }
}

/// <summary>
/// Simple texture renderer for FontStashSharp integration with OpenGL
/// </summary>
public class TextureRenderer : IFontStashRenderer, ITexture2DManager, IDisposable
{
    private int _shaderProgram;
    private int _vao;
    private int _vbo;
    private int _ebo;
    private int _projectionLoc;
    private int _textureLoc;

    private List<VertexPositionColorTexture> _vertices = new();
    private List<int> _indices = new();
    private Dictionary<object, int> _textureIds = new();

    private int _windowWidth;
    private int _windowHeight;
    private bool _isDrawing;

    public TextureRenderer()
    {
        CreateShaders();
        CreateBuffers();
    }

    private void CreateShaders()
    {
        string vertexSource = @"
            #version 330 core
            layout(location = 0) in vec2 aPosition;
            layout(location = 1) in vec4 aColor;
            layout(location = 2) in vec2 aTexCoord;

            uniform mat4 projection;

            out vec4 vColor;
            out vec2 vTexCoord;

            void main()
            {
                gl_Position = projection * vec4(aPosition, 0.0, 1.0);
                vColor = aColor;
                vTexCoord = aTexCoord;
            }
        ";

        string fragmentSource = @"
            #version 330 core
            in vec4 vColor;
            in vec2 vTexCoord;

            uniform sampler2D uTexture;

            out vec4 FragColor;

            void main()
            {
                vec4 texColor = texture(uTexture, vTexCoord);
                FragColor = vColor * texColor;
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

        _projectionLoc = GL.GetUniformLocation(_shaderProgram, "projection");
        _textureLoc = GL.GetUniformLocation(_shaderProgram, "uTexture");
    }

    private void CreateBuffers()
    {
        _vao = GL.GenVertexArray();
        _vbo = GL.GenBuffer();
        _ebo = GL.GenBuffer();

        GL.BindVertexArray(_vao);

        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);
        GL.BindBuffer(BufferTarget.ElementArrayBuffer, _ebo);

        int stride = 8 * sizeof(float);

        // Position
        GL.VertexAttribPointer(0, 2, VertexAttribPointerType.Float, false, stride, 0);
        GL.EnableVertexAttribArray(0);

        // Color
        GL.VertexAttribPointer(1, 4, VertexAttribPointerType.Float, false, stride, 2 * sizeof(float));
        GL.EnableVertexAttribArray(1);

        // TexCoord
        GL.VertexAttribPointer(2, 2, VertexAttribPointerType.Float, false, stride, 6 * sizeof(float));
        GL.EnableVertexAttribArray(2);

        GL.BindVertexArray(0);
    }

    public void Begin(int windowWidth, int windowHeight)
    {
        _windowWidth = windowWidth;
        _windowHeight = windowHeight;
        _vertices.Clear();
        _indices.Clear();
        _isDrawing = true;
    }

    public void End()
    {
        if (!_isDrawing || _vertices.Count == 0) return;

        GL.Enable(EnableCap.Blend);
        GL.BlendFunc(BlendingFactor.SrcAlpha, BlendingFactor.OneMinusSrcAlpha);

        GL.UseProgram(_shaderProgram);

        // Set orthographic projection for screen-space rendering
        Matrix4 projection = Matrix4.CreateOrthographicOffCenter(0, _windowWidth, _windowHeight, 0, -1, 1);
        GL.UniformMatrix4(_projectionLoc, false, ref projection);

        GL.ActiveTexture(TextureUnit.Texture0);
        GL.Uniform1(_textureLoc, 0);

        // Upload vertex data
        float[] vertexData = new float[_vertices.Count * 8];
        for (int i = 0; i < _vertices.Count; i++)
        {
            var v = _vertices[i];
            vertexData[i * 8 + 0] = v.Position.X;
            vertexData[i * 8 + 1] = v.Position.Y;
            vertexData[i * 8 + 2] = v.Color.R / 255f;
            vertexData[i * 8 + 3] = v.Color.G / 255f;
            vertexData[i * 8 + 4] = v.Color.B / 255f;
            vertexData[i * 8 + 5] = v.Color.A / 255f;
            vertexData[i * 8 + 6] = v.TextureCoordinate.X;
            vertexData[i * 8 + 7] = v.TextureCoordinate.Y;
        }

        GL.BindVertexArray(_vao);

        GL.BindBuffer(BufferTarget.ArrayBuffer, _vbo);
        GL.BufferData(BufferTarget.ArrayBuffer, vertexData.Length * sizeof(float), vertexData, BufferUsageHint.DynamicDraw);

        int[] indexData = _indices.ToArray();
        GL.BindBuffer(BufferTarget.ElementArrayBuffer, _ebo);
        GL.BufferData(BufferTarget.ElementArrayBuffer, indexData.Length * sizeof(int), indexData, BufferUsageHint.DynamicDraw);

        GL.DrawElements(PrimitiveType.Triangles, _indices.Count, DrawElementsType.UnsignedInt, 0);

        GL.BindVertexArray(0);
        GL.Disable(EnableCap.Blend);

        _isDrawing = false;
    }

    public ITexture2DManager TextureManager => this;

    public object CreateTexture(int width, int height)
    {
        return new Texture2D(width, height);
    }

    public System.Drawing.Point GetTextureSize(object texture)
    {
        if (texture is Texture2D tex)
        {
            return new System.Drawing.Point(tex.Width, tex.Height);
        }

        return System.Drawing.Point.Empty;
    }

    public void SetTextureData(object texture, System.Drawing.Rectangle bounds, byte[] data)
    {
        if (texture is Texture2D tex)
        {
            tex.SetData(bounds, data);
        }
    }

    public void Draw(object texture, System.Numerics.Vector2 pos, System.Drawing.Rectangle? src, FSColor color, float rotation,
        System.Numerics.Vector2 origin, float scale)
    {
        var scaleVector = new System.Numerics.Vector2(scale, scale);
        if (!_textureIds.TryGetValue(texture, out int textureId))
        {
            // Create OpenGL texture from FontStashSharp texture
            if (texture is Texture2D fsTexture)
            {
                textureId = CreateGlTexture(fsTexture);
                _textureIds[texture] = textureId;
            }
        }

        GL.BindTexture(TextureTarget.Texture2D, textureId);

        var srcRect = src ?? new System.Drawing.Rectangle(0, 0, 1, 1);

        float x = pos.X - origin.X * scaleVector.X;
        float y = pos.Y - origin.Y * scaleVector.Y;
        float w = srcRect.Width * scaleVector.X;
        float h = srcRect.Height * scaleVector.Y;

        // UV coordinates
        float u0 = 0, v0 = 0, u1 = 1, v1 = 1;
        if (texture is Texture2D tex)
        {
            u0 = (float)srcRect.X / tex.Width;
            v0 = (float)srcRect.Y / tex.Height;
            u1 = (float)(srcRect.X + srcRect.Width) / tex.Width;
            v1 = (float)(srcRect.Y + srcRect.Height) / tex.Height;
        }

        int startIndex = _vertices.Count;

        _vertices.Add(new VertexPositionColorTexture(new System.Numerics.Vector2(x, y), color, new System.Numerics.Vector2(u0, v0)));
        _vertices.Add(new VertexPositionColorTexture(new System.Numerics.Vector2(x + w, y), color, new System.Numerics.Vector2(u1, v0)));
        _vertices.Add(new VertexPositionColorTexture(new System.Numerics.Vector2(x + w, y + h), color, new System.Numerics.Vector2(u1, v1)));
        _vertices.Add(new VertexPositionColorTexture(new System.Numerics.Vector2(x, y + h), color, new System.Numerics.Vector2(u0, v1)));

        _indices.Add(startIndex);
        _indices.Add(startIndex + 1);
        _indices.Add(startIndex + 2);
        _indices.Add(startIndex);
        _indices.Add(startIndex + 2);
        _indices.Add(startIndex + 3);
    }

    private int CreateGlTexture(Texture2D fsTexture)
    {
        int textureId = GL.GenTexture();
        GL.BindTexture(TextureTarget.Texture2D, textureId);

        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMinFilter, (int)TextureMinFilter.Linear);
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);

        // Get pixel data from FontStashSharp texture
        byte[] data = new byte[fsTexture.Width * fsTexture.Height * 4];
        fsTexture.GetData(data);

        GL.TexImage2D(TextureTarget.Texture2D, 0, PixelInternalFormat.Rgba,
            fsTexture.Width, fsTexture.Height, 0,
            PixelFormat.Rgba, PixelType.UnsignedByte, data);

        return textureId;
    }

    public void Dispose()
    {
        foreach (var textureId in _textureIds.Values)
        {
            GL.DeleteTexture(textureId);
        }

        GL.DeleteVertexArray(_vao);
        GL.DeleteBuffer(_vbo);
        GL.DeleteBuffer(_ebo);
        GL.DeleteProgram(_shaderProgram);
    }
}

public struct VertexPositionColorTexture
{
    public System.Numerics.Vector2 Position;
    public FSColor Color;
    public System.Numerics.Vector2 TextureCoordinate;

    public VertexPositionColorTexture(System.Numerics.Vector2 position, FSColor color, System.Numerics.Vector2 texCoord)
    {
        Position = position;
        Color = color;
        TextureCoordinate = texCoord;
    }
}

/// <summary>
/// FontStashSharp texture wrapper
/// </summary>
public class Texture2D
{
    private int _textureId;
    private byte[]? _pixelData;

    public int Width { get; }
    public int Height { get; }

    public Texture2D(int width, int height)
    {
        Width = width;
        Height = height;
        _pixelData = new byte[width * height * 4];

        _textureId = GL.GenTexture();
        GL.BindTexture(TextureTarget.Texture2D, _textureId);
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMinFilter, (int)TextureMinFilter.Linear);
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);
        GL.TexImage2D(TextureTarget.Texture2D, 0, PixelInternalFormat.Rgba,
            width, height, 0, PixelFormat.Rgba, PixelType.UnsignedByte, _pixelData);
    }

    public void SetData(System.Drawing.Rectangle bounds, byte[] data)
    {
        if (_pixelData == null) return;

        // Copy data to pixel buffer
        for (int y = 0; y < bounds.Height; y++)
        {
            for (int x = 0; x < bounds.Width; x++)
            {
                int srcIdx = (y * bounds.Width + x) * 4;
                int dstIdx = ((bounds.Y + y) * Width + (bounds.X + x)) * 4;

                if (srcIdx + 3 < data.Length && dstIdx + 3 < _pixelData.Length)
                {
                    _pixelData[dstIdx] = data[srcIdx];
                    _pixelData[dstIdx + 1] = data[srcIdx + 1];
                    _pixelData[dstIdx + 2] = data[srcIdx + 2];
                    _pixelData[dstIdx + 3] = data[srcIdx + 3];
                }
            }
        }

        // Update OpenGL texture
        GL.BindTexture(TextureTarget.Texture2D, _textureId);
        GL.TexSubImage2D(TextureTarget.Texture2D, 0, bounds.X, bounds.Y, bounds.Width, bounds.Height,
            PixelFormat.Rgba, PixelType.UnsignedByte, data);
    }

    public void GetData(byte[] data)
    {
        if (_pixelData != null)
        {
            Array.Copy(_pixelData, data, Math.Min(data.Length, _pixelData.Length));
        }
    }

    public void Dispose()
    {
        GL.DeleteTexture(_textureId);
        _pixelData = null;
    }
}
