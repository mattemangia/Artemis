using System;
using System.Collections.Generic;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Spatial hash grid for broad-phase collision detection optimization.
    /// Divides space into cells and only checks collisions within same/nearby cells.
    /// </summary>
    public class SpatialHashGrid2D
    {
        private Dictionary<(int, int), List<RigidBody2D>> _grid;
        private double _cellSize;
        private double _invCellSize;

        public SpatialHashGrid2D(double cellSize = 10.0)
        {
            _cellSize = cellSize;
            _invCellSize = 1.0 / cellSize;
            _grid = new Dictionary<(int, int), List<RigidBody2D>>();
        }

        public void Clear()
        {
            foreach (var cell in _grid.Values)
            {
                cell.Clear();
            }
            _grid.Clear();
        }

        public void Insert(RigidBody2D body)
        {
            var cells = GetCells(body);
            foreach (var cell in cells)
            {
                if (!_grid.ContainsKey(cell))
                {
                    _grid[cell] = new List<RigidBody2D>();
                }
                _grid[cell].Add(body);
            }
        }

        public List<RigidBody2D> QueryNearby(RigidBody2D body)
        {
            var nearby = new HashSet<RigidBody2D>();
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

            return new List<RigidBody2D>(nearby);
        }

        public List<RigidBody2D> QueryRegion(Vector2D min, Vector2D max)
        {
            var bodies = new HashSet<RigidBody2D>();

            int minCellX = (int)Math.Floor(min.X * _invCellSize);
            int minCellY = (int)Math.Floor(min.Y * _invCellSize);
            int maxCellX = (int)Math.Floor(max.X * _invCellSize);
            int maxCellY = (int)Math.Floor(max.Y * _invCellSize);

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

            return new List<RigidBody2D>(bodies);
        }

        private List<(int, int)> GetCells(RigidBody2D body)
        {
            var cells = new List<(int, int)>();
            var aabb = body.AABB;

            int minCellX = (int)Math.Floor(aabb.Min.X * _invCellSize);
            int minCellY = (int)Math.Floor(aabb.Min.Y * _invCellSize);
            int maxCellX = (int)Math.Floor(aabb.Max.X * _invCellSize);
            int maxCellY = (int)Math.Floor(aabb.Max.Y * _invCellSize);

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
    /// Quadtree for spatial partitioning - alternative to hash grid.
    /// Better for non-uniform distributions of objects.
    /// </summary>
    public class QuadTree2D
    {
        private const int MaxObjectsPerNode = 4;
        private const int MaxDepth = 8;

        private int _depth;
        private List<RigidBody2D> _bodies;
        private AABB2D _bounds;
        private QuadTree2D?[] _children;
        private bool _isDivided;

        public QuadTree2D(AABB2D bounds, int depth = 0)
        {
            _bounds = bounds;
            _depth = depth;
            _bodies = new List<RigidBody2D>();
            _children = new QuadTree2D[4];
            _isDivided = false;
        }

        public void Clear()
        {
            _bodies.Clear();
            if (_isDivided)
            {
                for (int i = 0; i < 4; i++)
                {
                    _children[i]?.Clear();
                    _children[i] = null;
                }
                _isDivided = false;
            }
        }

        public void Insert(RigidBody2D body)
        {
            if (!_bounds.Contains(body.Position))
                return;

            if (_isDivided)
            {
                int index = GetChildIndex(body.Position);
                if (index != -1 && _children[index] != null)
                {
                    _children[index]!.Insert(body);
                    return;
                }
            }

            _bodies.Add(body);

            if (_bodies.Count > MaxObjectsPerNode && _depth < MaxDepth && !_isDivided)
            {
                Subdivide();

                int i = 0;
                while (i < _bodies.Count)
                {
                    int index = GetChildIndex(_bodies[i].Position);
                    if (index != -1 && _children[index] != null)
                    {
                        _children[index]!.Insert(_bodies[i]);
                        _bodies.RemoveAt(i);
                    }
                    else
                    {
                        i++;
                    }
                }
            }
        }

        public List<RigidBody2D> Query(AABB2D range)
        {
            var result = new List<RigidBody2D>();

            if (!_bounds.Overlaps(range))
                return result;

            foreach (var body in _bodies)
            {
                if (range.Contains(body.Position))
                {
                    result.Add(body);
                }
            }

            if (_isDivided)
            {
                for (int i = 0; i < 4; i++)
                {
                    if (_children[i] != null)
                    {
                        result.AddRange(_children[i]!.Query(range));
                    }
                }
            }

            return result;
        }

        public List<RigidBody2D> QueryRadius(Vector2D center, double radius)
        {
            // Create AABB for the circle
            var aabb = new AABB2D(
                center - new Vector2D(radius, radius),
                center + new Vector2D(radius, radius)
            );

            var candidates = Query(aabb);
            var result = new List<RigidBody2D>();

            double radiusSq = radius * radius;
            foreach (var body in candidates)
            {
                if ((body.Position - center).MagnitudeSquared <= radiusSq)
                {
                    result.Add(body);
                }
            }

            return result;
        }

        private void Subdivide()
        {
            double halfW = _bounds.Width / 2;
            double halfH = _bounds.Height / 2;
            Vector2D center = _bounds.Center;

            _children[0] = new QuadTree2D(new AABB2D(
                new Vector2D(center.X, center.Y),
                new Vector2D(center.X + halfW, center.Y + halfH)), _depth + 1);

            _children[1] = new QuadTree2D(new AABB2D(
                new Vector2D(center.X - halfW, center.Y),
                new Vector2D(center.X, center.Y + halfH)), _depth + 1);

            _children[2] = new QuadTree2D(new AABB2D(
                new Vector2D(center.X - halfW, center.Y - halfH),
                new Vector2D(center.X, center.Y)), _depth + 1);

            _children[3] = new QuadTree2D(new AABB2D(
                new Vector2D(center.X, center.Y - halfH),
                new Vector2D(center.X + halfW, center.Y)), _depth + 1);

            _isDivided = true;
        }

        private int GetChildIndex(Vector2D position)
        {
            Vector2D center = _bounds.Center;

            bool right = position.X >= center.X;
            bool top = position.Y >= center.Y;

            if (right && top) return 0;
            if (!right && top) return 1;
            if (!right && !top) return 2;
            if (right && !top) return 3;

            return -1;
        }
    }
}
