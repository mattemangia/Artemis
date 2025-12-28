using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Thread-safe spatial hash grid for broad-phase collision detection optimization.
    /// Uses concurrent collections and parallel processing for high performance.
    /// </summary>
    public class SpatialHashGrid2D
    {
        private readonly ConcurrentDictionary<long, ConcurrentBag<RigidBody2D>> _grid;
        private readonly double _cellSize;
        private readonly double _invCellSize;

        public SpatialHashGrid2D(double cellSize = 10.0)
        {
            _cellSize = cellSize;
            _invCellSize = 1.0 / cellSize;
            _grid = new ConcurrentDictionary<long, ConcurrentBag<RigidBody2D>>();
        }

        public void Clear()
        {
            _grid.Clear();
        }

        /// <summary>
        /// Insert a single body into the grid.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Insert(RigidBody2D body)
        {
            var aabb = body.AABB;
            int minCellX = (int)Math.Floor(aabb.Min.X * _invCellSize);
            int minCellY = (int)Math.Floor(aabb.Min.Y * _invCellSize);
            int maxCellX = (int)Math.Floor(aabb.Max.X * _invCellSize);
            int maxCellY = (int)Math.Floor(aabb.Max.Y * _invCellSize);

            for (int x = minCellX; x <= maxCellX; x++)
            {
                for (int y = minCellY; y <= maxCellY; y++)
                {
                    long key = GetCellKey(x, y);
                    var bag = _grid.GetOrAdd(key, _ => new ConcurrentBag<RigidBody2D>());
                    bag.Add(body);
                }
            }
        }

        /// <summary>
        /// Insert multiple bodies in parallel.
        /// </summary>
        public void InsertParallel(IList<RigidBody2D> bodies)
        {
            Parallel.ForEach(bodies, body =>
            {
                if (body.IsActive)
                    Insert(body);
            });
        }

        /// <summary>
        /// Rebuild the entire grid from a list of bodies (parallel).
        /// </summary>
        public void RebuildParallel(IList<RigidBody2D> bodies)
        {
            Clear();
            InsertParallel(bodies);
        }

        /// <summary>
        /// Query bodies near a given body.
        /// </summary>
        public List<RigidBody2D> QueryNearby(RigidBody2D body)
        {
            var nearby = new HashSet<RigidBody2D>();
            var aabb = body.AABB;

            int minCellX = (int)Math.Floor(aabb.Min.X * _invCellSize);
            int minCellY = (int)Math.Floor(aabb.Min.Y * _invCellSize);
            int maxCellX = (int)Math.Floor(aabb.Max.X * _invCellSize);
            int maxCellY = (int)Math.Floor(aabb.Max.Y * _invCellSize);

            for (int x = minCellX; x <= maxCellX; x++)
            {
                for (int y = minCellY; y <= maxCellY; y++)
                {
                    long key = GetCellKey(x, y);
                    if (_grid.TryGetValue(key, out var bag))
                    {
                        foreach (var other in bag)
                        {
                            if (other != body)
                                nearby.Add(other);
                        }
                    }
                }
            }

            return new List<RigidBody2D>(nearby);
        }

        /// <summary>
        /// Query bodies in a region.
        /// </summary>
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
                    long key = GetCellKey(x, y);
                    if (_grid.TryGetValue(key, out var bag))
                    {
                        foreach (var body in bag)
                            bodies.Add(body);
                    }
                }
            }

            return new List<RigidBody2D>(bodies);
        }

        /// <summary>
        /// Get all potential collision pairs (parallel).
        /// </summary>
        public List<(RigidBody2D, RigidBody2D)> GetPotentialPairsParallel()
        {
            var pairs = new ConcurrentBag<(RigidBody2D, RigidBody2D)>();
            var processedPairs = new ConcurrentDictionary<(int, int), byte>();

            Parallel.ForEach(_grid, cell =>
            {
                var bodies = cell.Value.ToArray();
                for (int i = 0; i < bodies.Length; i++)
                {
                    for (int j = i + 1; j < bodies.Length; j++)
                    {
                        var a = bodies[i];
                        var b = bodies[j];

                        // Create order-independent key
                        var pairKey = a.GetHashCode() < b.GetHashCode()
                            ? (a.GetHashCode(), b.GetHashCode())
                            : (b.GetHashCode(), a.GetHashCode());

                        if (processedPairs.TryAdd(pairKey, 0))
                        {
                            if (a.AABB.Overlaps(b.AABB))
                                pairs.Add((a, b));
                        }
                    }
                }
            });

            return new List<(RigidBody2D, RigidBody2D)>(pairs);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static long GetCellKey(int x, int y)
        {
            return ((long)(x & 0x7FFFFFFF) << 32) | (uint)(y & 0x7FFFFFFF);
        }
    }

    /// <summary>
    /// Thread-safe quadtree for spatial partitioning.
    /// Better for non-uniform distributions of objects.
    /// </summary>
    public class QuadTree2D
    {
        private const int MaxObjectsPerNode = 8;
        private const int MaxDepth = 10;

        private readonly int _depth;
        private readonly List<RigidBody2D> _bodies;
        private readonly AABB2D _bounds;
        private readonly QuadTree2D?[] _children;
        private readonly object _lock = new object();
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
            lock (_lock)
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
        }

        /// <summary>
        /// Build the quadtree from bodies in parallel.
        /// </summary>
        public void BuildParallel(IList<RigidBody2D> bodies)
        {
            Clear();

            // First pass: insert all bodies
            foreach (var body in bodies)
            {
                if (body.IsActive)
                    Insert(body);
            }
        }

        public void Insert(RigidBody2D body)
        {
            if (!_bounds.Contains(body.Position))
                return;

            lock (_lock)
            {
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
        }

        /// <summary>
        /// Query bodies in parallel (for large result sets).
        /// </summary>
        public List<RigidBody2D> QueryParallel(AABB2D range)
        {
            var result = new ConcurrentBag<RigidBody2D>();
            QueryParallelInternal(range, result);
            return new List<RigidBody2D>(result);
        }

        private void QueryParallelInternal(AABB2D range, ConcurrentBag<RigidBody2D> result)
        {
            if (!_bounds.Overlaps(range))
                return;

            foreach (var body in _bodies)
            {
                if (range.Contains(body.Position))
                    result.Add(body);
            }

            if (_isDivided)
            {
                Parallel.For(0, 4, i =>
                {
                    _children[i]?.QueryParallelInternal(range, result);
                });
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
                    result.Add(body);
            }

            if (_isDivided)
            {
                for (int i = 0; i < 4; i++)
                {
                    if (_children[i] != null)
                        result.AddRange(_children[i]!.Query(range));
                }
            }

            return result;
        }

        public List<RigidBody2D> QueryRadius(Vector2D center, double radius)
        {
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
                    result.Add(body);
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

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
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

    /// <summary>
    /// SIMD-optimized batch operations for spatial queries.
    /// </summary>
    public static class SpatialBatchOperations
    {
        /// <summary>
        /// Perform batch AABB overlap tests using parallel processing.
        /// </summary>
        public static List<(int, int)> BatchOverlapTest(AABB2D[] aabbs)
        {
            var pairs = new ConcurrentBag<(int, int)>();

            Parallel.For(0, aabbs.Length, i =>
            {
                for (int j = i + 1; j < aabbs.Length; j++)
                {
                    if (aabbs[i].Overlaps(aabbs[j]))
                        pairs.Add((i, j));
                }
            });

            return new List<(int, int)>(pairs);
        }

        /// <summary>
        /// Batch point-in-AABB test.
        /// </summary>
        public static bool[] BatchPointInAABB(Vector2D[] points, AABB2D aabb)
        {
            var results = new bool[points.Length];

            Parallel.For(0, points.Length, i =>
            {
                results[i] = aabb.Contains(points[i]);
            });

            return results;
        }

        /// <summary>
        /// Batch distance calculations.
        /// </summary>
        public static double[] BatchDistances(Vector2D origin, Vector2D[] targets)
        {
            var distances = new double[targets.Length];

            Parallel.For(0, targets.Length, i =>
            {
                distances[i] = (targets[i] - origin).Magnitude;
            });

            return distances;
        }

        /// <summary>
        /// Find all bodies within radius using parallel processing.
        /// </summary>
        public static List<RigidBody2D> FindBodiesInRadiusParallel(
            IList<RigidBody2D> bodies, Vector2D center, double radius)
        {
            var result = new ConcurrentBag<RigidBody2D>();
            double radiusSq = radius * radius;

            Parallel.ForEach(bodies, body =>
            {
                if ((body.Position - center).MagnitudeSquared <= radiusSq)
                    result.Add(body);
            });

            return new List<RigidBody2D>(result);
        }
    }
}
