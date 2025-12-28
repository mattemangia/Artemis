
namespace PhysicsCatapultDemo;

public class Renderer
{
    private int _width;
    private int _height;
    private char[,] _buffer;
    private ConsoleColor[,] _colorBuffer;
    private float _scale;
    private Vector2 _offset;

    public Renderer(int width, int height, float scale = 2f)
    {
        _width = width;
        _height = height;
        _scale = scale;
        _offset = new Vector2(5, 5);
        _buffer = new char[height, width];
        _colorBuffer = new ConsoleColor[height, width];
    }

    public void Clear()
    {
        for (int y = 0; y < _height; y++)
        {
            for (int x = 0; x < _width; x++)
            {
                _buffer[y, x] = ' ';
                _colorBuffer[y, x] = ConsoleColor.Black;
            }
        }
    }

    public void DrawLine(Vector2 start, Vector2 end, char c, ConsoleColor color)
    {
        Vector2 startScreen = WorldToScreen(start);
        Vector2 endScreen = WorldToScreen(end);

        int x0 = (int)startScreen.X;
        int y0 = (int)startScreen.Y;
        int x1 = (int)endScreen.X;
        int y1 = (int)endScreen.Y;

        int dx = Math.Abs(x1 - x0);
        int dy = Math.Abs(y1 - y0);
        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;
        int err = dx - dy;

        while (true)
        {
            SetPixel(x0, y0, c, color);

            if (x0 == x1 && y0 == y1) break;

            int e2 = 2 * err;
            if (e2 > -dy)
            {
                err -= dy;
                x0 += sx;
            }
            if (e2 < dx)
            {
                err += dx;
                y0 += sy;
            }
        }
    }

    public void DrawCircle(Vector2 center, double radius, char c, ConsoleColor color)
    {
        Vector2 screenPos = WorldToScreen(center);
        int screenRadius = (int)(radius * _scale);

        for (int dy = -screenRadius; dy <= screenRadius; dy++)
        {
            for (int dx = -screenRadius; dx <= screenRadius; dx++)
            {
                if (dx * dx + dy * dy <= screenRadius * screenRadius)
                {
                    int x = (int)screenPos.X + dx;
                    int y = (int)screenPos.Y + dy;
                    SetPixel(x, y, c, color);
                }
            }
        }
    }

    public void DrawBox(Vector2 center, double width, double height, char c, ConsoleColor color)
    {
        Vector2 screenPos = WorldToScreen(center);
        int screenWidth = (int)(width * _scale);
        int screenHeight = (int)(height * _scale);

        int startX = (int)screenPos.X - screenWidth / 2;
        int startY = (int)screenPos.Y - screenHeight / 2;

        for (int dy = 0; dy < screenHeight; dy++)
        {
            for (int dx = 0; dx < screenWidth; dx++)
            {
                SetPixel(startX + dx, startY + dy, c, color);
            }
        }
    }

    public void DrawText(int x, int y, string text, ConsoleColor color = ConsoleColor.White)
    {
        for (int i = 0; i < text.Length && x + i < _width; i++)
        {
            SetPixel(x + i, y, text[i], color);
        }
    }

    private void SetPixel(int x, int y, char c, ConsoleColor color)
    {
        if (x >= 0 && x < _width && y >= 0 && y < _height)
        {
            _buffer[y, x] = c;
            _colorBuffer[y, x] = color;
        }
    }

    private Vector2 WorldToScreen(Vector2 worldPos)
    {
        return new Vector2(
            (worldPos.X + _offset.X) * _scale,
            _height - (worldPos.Y + _offset.Y) * _scale
        );
    }

    public void Present()
    {
        Console.SetCursorPosition(0, 0);

        for (int y = 0; y < _height; y++)
        {
            Console.SetCursorPosition(0, y);
            for (int x = 0; x < _width; x++)
            {
                Console.ForegroundColor = _colorBuffer[y, x];
                Console.Write(_buffer[y, x]);
            }
            // Don't use WriteLine - use SetCursorPosition to avoid scrolling
        }

        Console.ResetColor();
    }
}
