namespace ArtemisEngine;

/// <summary>
/// Spatial hash grid for broad-phase collision detection optimization
/// Divides space into cells and only checks collisions within same/nearby cells
/// </summary>
public class SpatialHashGrid
{
    private Dictionary<(int, int), List<RigidBody>> _grid;
    private float _cellSize;

    public SpatialHashGrid(float cellSize = 10f)
    {
        _cellSize = cellSize;
        _grid = new Dictionary<(int, int), List<RigidBody>>();
    }

    public void Clear()
    {
        foreach (var cell in _grid.Values)
        {
            cell.Clear();
        }
        _grid.Clear();
    }

    public void Insert(RigidBody body)
    {
        var cells = GetCells(body);
        foreach (var cell in cells)
        {
            if (!_grid.ContainsKey(cell))
            {
                _grid[cell] = new List<RigidBody>();
            }
            _grid[cell].Add(body);
        }
    }

    public List<RigidBody> QueryNearby(RigidBody body)
    {
        var nearby = new HashSet<RigidBody>();
        var cells = GetCells(body);

        foreach (var cell in cells)
        {
            if (_grid.TryGetValue(cell, out var bodies))
            {
                foreach (var other in bodies)
                {
                    if (other != body)
                    {
                        nearby.Add(other);
                    }
                }
            }
        }

        return nearby.ToList();
    }

    public List<RigidBody> QueryRegion(Vector2 min, Vector2 max)
    {
        var bodies = new HashSet<RigidBody>();

        int minCellX = (int)MathF.Floor(min.X / _cellSize);
        int minCellY = (int)MathF.Floor(min.Y / _cellSize);
        int maxCellX = (int)MathF.Floor(max.X / _cellSize);
        int maxCellY = (int)MathF.Floor(max.Y / _cellSize);

        for (int x = minCellX; x <= maxCellX; x++)
        {
            for (int y = minCellY; y <= maxCellY; y++)
            {
                if (_grid.TryGetValue((x, y), out var cellBodies))
                {
                    foreach (var body in cellBodies)
                    {
                        bodies.Add(body);
                    }
                }
            }
        }

        return bodies.ToList();
    }

    private List<(int, int)> GetCells(RigidBody body)
    {
        var cells = new List<(int, int)>();
        Vector2 min, max;

        if (body.Shape is CircleShape circle)
        {
            min = new Vector2(body.Position.X - circle.Radius, body.Position.Y - circle.Radius);
            max = new Vector2(body.Position.X + circle.Radius, body.Position.Y + circle.Radius);
        }
        else if (body.Shape is BoxShape box)
        {
            float halfW = box.Width / 2;
            float halfH = box.Height / 2;
            min = new Vector2(body.Position.X - halfW, body.Position.Y - halfH);
            max = new Vector2(body.Position.X + halfW, body.Position.Y + halfH);
        }
        else
        {
            return cells;
        }

        int minCellX = (int)MathF.Floor(min.X / _cellSize);
        int minCellY = (int)MathF.Floor(min.Y / _cellSize);
        int maxCellX = (int)MathF.Floor(max.X / _cellSize);
        int maxCellY = (int)MathF.Floor(max.Y / _cellSize);

        for (int x = minCellX; x <= maxCellX; x++)
        {
            for (int y = minCellY; y <= maxCellY; y++)
            {
                cells.Add((x, y));
            }
        }

        return cells;
    }
}

/// <summary>
/// Quadtree for spatial partitioning - alternative to hash grid
/// Better for non-uniform distributions of objects
/// </summary>
public class QuadTree
{
    private const int MaxObjectsPerNode = 4;
    private const int MaxDepth = 5;

    private int _depth;
    private List<RigidBody> _bodies;
    private AABB _bounds;
    private QuadTree[] _children;

    public QuadTree(AABB bounds, int depth = 0)
    {
        _bounds = bounds;
        _depth = depth;
        _bodies = new List<RigidBody>();
        _children = new QuadTree[4];
    }

    public void Clear()
    {
        _bodies.Clear();
        for (int i = 0; i < 4; i++)
        {
            _children[i]?.Clear();
            _children[i] = null!;
        }
    }

    public void Insert(RigidBody body)
    {
        if (!_bounds.Contains(body.Position))
            return;

        if (_children[0] != null)
        {
            int index = GetChildIndex(body.Position);
            if (index != -1)
            {
                _children[index].Insert(body);
                return;
            }
        }

        _bodies.Add(body);

        if (_bodies.Count > MaxObjectsPerNode && _depth < MaxDepth)
        {
            if (_children[0] == null)
            {
                Subdivide();
            }

            int i = 0;
            while (i < _bodies.Count)
            {
                int index = GetChildIndex(_bodies[i].Position);
                if (index != -1)
                {
                    _children[index].Insert(_bodies[i]);
                    _bodies.RemoveAt(i);
                }
                else
                {
                    i++;
                }
            }
        }
    }

    public List<RigidBody> Query(AABB range)
    {
        var result = new List<RigidBody>();

        if (!_bounds.Intersects(range))
            return result;

        foreach (var body in _bodies)
        {
            if (range.Contains(body.Position))
            {
                result.Add(body);
            }
        }

        if (_children[0] != null)
        {
            for (int i = 0; i < 4; i++)
            {
                result.AddRange(_children[i].Query(range));
            }
        }

        return result;
    }

    private void Subdivide()
    {
        float halfW = _bounds.Width / 2;
        float halfH = _bounds.Height / 2;
        float x = _bounds.Center.X;
        float y = _bounds.Center.Y;

        _children[0] = new QuadTree(new AABB(new Vector2(x + halfW / 2, y + halfH / 2), halfW, halfH), _depth + 1);
        _children[1] = new QuadTree(new AABB(new Vector2(x - halfW / 2, y + halfH / 2), halfW, halfH), _depth + 1);
        _children[2] = new QuadTree(new AABB(new Vector2(x - halfW / 2, y - halfH / 2), halfW, halfH), _depth + 1);
        _children[3] = new QuadTree(new AABB(new Vector2(x + halfW / 2, y - halfH / 2), halfW, halfH), _depth + 1);
    }

    private int GetChildIndex(Vector2 position)
    {
        float midX = _bounds.Center.X;
        float midY = _bounds.Center.Y;

        bool right = position.X >= midX;
        bool top = position.Y >= midY;

        if (right && top) return 0;
        if (!right && top) return 1;
        if (!right && !top) return 2;
        if (right && !top) return 3;

        return -1;
    }
}

public struct AABB
{
    public Vector2 Center { get; set; }
    public float Width { get; set; }
    public float Height { get; set; }

    public AABB(Vector2 center, float width, float height)
    {
        Center = center;
        Width = width;
        Height = height;
    }

    public Vector2 Min => new Vector2(Center.X - Width / 2, Center.Y - Height / 2);
    public Vector2 Max => new Vector2(Center.X + Width / 2, Center.Y + Height / 2);

    public bool Contains(Vector2 point)
    {
        return point.X >= Min.X && point.X <= Max.X &&
               point.Y >= Min.Y && point.Y <= Max.Y;
    }

    public bool Intersects(AABB other)
    {
        return !(Min.X > other.Max.X || Max.X < other.Min.X ||
                 Min.Y > other.Max.Y || Max.Y < other.Min.Y);
    }
}
