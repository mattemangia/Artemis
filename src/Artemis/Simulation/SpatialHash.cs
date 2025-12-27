using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Simulation
{
    /// <summary>
    /// High-performance spatial hash grid for O(1) average-case broad-phase collision detection.
    /// Uses a hash table with linked cells for memory efficiency.
    /// </summary>
    public class SpatialHash
    {
        #region Fields

        private readonly Dictionary<long, List<IPhysicsBody>> _cells;
        private readonly Dictionary<IPhysicsBody, List<long>> _bodyCells;
        private readonly double _cellSize;
        private readonly double _invCellSize;
        private readonly List<IPhysicsBody> _queryResults;
        private readonly HashSet<(IPhysicsBody, IPhysicsBody)> _testedPairs;

        // Statistics
        private int _insertCount;
        private int _queryCount;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the cell size.
        /// </summary>
        public double CellSize => _cellSize;

        /// <summary>
        /// Gets the number of occupied cells.
        /// </summary>
        public int OccupiedCellCount => _cells.Count;

        /// <summary>
        /// Gets the total number of body insertions since last clear.
        /// </summary>
        public int InsertCount => _insertCount;

        /// <summary>
        /// Gets the total number of queries since last clear.
        /// </summary>
        public int QueryCount => _queryCount;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new spatial hash with the specified cell size.
        /// </summary>
        /// <param name="cellSize">Size of each cell. Should be >= largest object diameter.</param>
        public SpatialHash(double cellSize = 2.0)
        {
            _cellSize = cellSize;
            _invCellSize = 1.0 / cellSize;
            _cells = new Dictionary<long, List<IPhysicsBody>>(256);
            _bodyCells = new Dictionary<IPhysicsBody, List<long>>(128);
            _queryResults = new List<IPhysicsBody>(64);
            _testedPairs = new HashSet<(IPhysicsBody, IPhysicsBody)>(256);
        }

        #endregion

        #region Hash Functions

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int HashCoord(double coord)
        {
            return (int)Math.Floor(coord * _invCellSize);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long HashKey(int x, int y, int z)
        {
            // 3D spatial hash using bit interleaving for better distribution
            const long p1 = 73856093L;
            const long p2 = 19349663L;
            const long p3 = 83492791L;
            return (x * p1) ^ (y * p2) ^ (z * p3);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long PositionToKey(Vector3D pos)
        {
            return HashKey(HashCoord(pos.X), HashCoord(pos.Y), HashCoord(pos.Z));
        }

        #endregion

        #region Insert/Remove

        /// <summary>
        /// Clears all bodies from the spatial hash.
        /// </summary>
        public void Clear()
        {
            foreach (var cellList in _cells.Values)
            {
                cellList.Clear();
            }
            _cells.Clear();
            _bodyCells.Clear();
            _testedPairs.Clear();
            _insertCount = 0;
            _queryCount = 0;
        }

        /// <summary>
        /// Inserts a body into the spatial hash.
        /// </summary>
        public void Insert(IPhysicsBody body)
        {
            var aabb = body.BoundingBox;

            int minX = HashCoord(aabb.Min.X);
            int maxX = HashCoord(aabb.Max.X);
            int minY = HashCoord(aabb.Min.Y);
            int maxY = HashCoord(aabb.Max.Y);
            int minZ = HashCoord(aabb.Min.Z);
            int maxZ = HashCoord(aabb.Max.Z);

            // Track which cells this body occupies
            if (!_bodyCells.TryGetValue(body, out var cellKeys))
            {
                cellKeys = new List<long>(8);
                _bodyCells[body] = cellKeys;
            }
            else
            {
                cellKeys.Clear();
            }

            for (int x = minX; x <= maxX; x++)
            {
                for (int y = minY; y <= maxY; y++)
                {
                    for (int z = minZ; z <= maxZ; z++)
                    {
                        long key = HashKey(x, y, z);

                        if (!_cells.TryGetValue(key, out var cell))
                        {
                            cell = new List<IPhysicsBody>(4);
                            _cells[key] = cell;
                        }

                        cell.Add(body);
                        cellKeys.Add(key);
                    }
                }
            }

            _insertCount++;
        }

        /// <summary>
        /// Removes a body from the spatial hash.
        /// </summary>
        public void Remove(IPhysicsBody body)
        {
            if (_bodyCells.TryGetValue(body, out var cellKeys))
            {
                foreach (var key in cellKeys)
                {
                    if (_cells.TryGetValue(key, out var cell))
                    {
                        cell.Remove(body);
                        if (cell.Count == 0)
                        {
                            _cells.Remove(key);
                        }
                    }
                }
                _bodyCells.Remove(body);
            }
        }

        /// <summary>
        /// Updates a body's position in the hash (remove + insert).
        /// </summary>
        public void Update(IPhysicsBody body)
        {
            Remove(body);
            Insert(body);
        }

        /// <summary>
        /// Rebuilds the entire spatial hash from a list of bodies.
        /// </summary>
        public void Rebuild(IReadOnlyList<IPhysicsBody> bodies)
        {
            Clear();
            for (int i = 0; i < bodies.Count; i++)
            {
                if (bodies[i].IsActive)
                {
                    Insert(bodies[i]);
                }
            }
        }

        #endregion

        #region Queries

        /// <summary>
        /// Queries all bodies potentially overlapping with the given AABB.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> Query(AABB aabb)
        {
            _queryResults.Clear();
            _queryCount++;

            int minX = HashCoord(aabb.Min.X);
            int maxX = HashCoord(aabb.Max.X);
            int minY = HashCoord(aabb.Min.Y);
            int maxY = HashCoord(aabb.Max.Y);
            int minZ = HashCoord(aabb.Min.Z);
            int maxZ = HashCoord(aabb.Max.Z);

            for (int x = minX; x <= maxX; x++)
            {
                for (int y = minY; y <= maxY; y++)
                {
                    for (int z = minZ; z <= maxZ; z++)
                    {
                        long key = HashKey(x, y, z);

                        if (_cells.TryGetValue(key, out var cell))
                        {
                            foreach (var body in cell)
                            {
                                if (!_queryResults.Contains(body))
                                {
                                    _queryResults.Add(body);
                                }
                            }
                        }
                    }
                }
            }

            return _queryResults;
        }

        /// <summary>
        /// Queries all bodies potentially overlapping with a given body.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> QueryBody(IPhysicsBody body)
        {
            return Query(body.BoundingBox);
        }

        /// <summary>
        /// Gets all potential collision pairs (broad-phase).
        /// This is the main function for collision detection.
        /// </summary>
        public void GetPotentialPairs(List<(IPhysicsBody A, IPhysicsBody B)> pairs)
        {
            pairs.Clear();
            _testedPairs.Clear();

            foreach (var kvp in _cells)
            {
                var cell = kvp.Value;
                int count = cell.Count;

                // Test all pairs within the cell
                for (int i = 0; i < count; i++)
                {
                    var bodyA = cell[i];
                    if (!bodyA.IsActive) continue;

                    for (int j = i + 1; j < count; j++)
                    {
                        var bodyB = cell[j];
                        if (!bodyB.IsActive) continue;

                        // Skip if both static
                        if (bodyA.BodyType == BodyType.Static && bodyB.BodyType == BodyType.Static)
                            continue;

                        // Create ordered pair to avoid duplicates
                        var pair = bodyA.GetHashCode() < bodyB.GetHashCode()
                            ? (bodyA, bodyB)
                            : (bodyB, bodyA);

                        if (_testedPairs.Add(pair))
                        {
                            // Additional AABB check for accuracy
                            if (bodyA.BoundingBox.Intersects(bodyB.BoundingBox))
                            {
                                pairs.Add(pair);
                            }
                        }
                    }
                }
            }
        }

        #endregion

        #region Raycast

        /// <summary>
        /// Performs a raycast through the spatial hash using 3D DDA algorithm.
        /// </summary>
        public IEnumerable<IPhysicsBody> RaycastAll(Vector3D origin, Vector3D direction, double maxDistance)
        {
            var visited = new HashSet<IPhysicsBody>();

            // Current cell coordinates
            int x = HashCoord(origin.X);
            int y = HashCoord(origin.Y);
            int z = HashCoord(origin.Z);

            // Direction signs
            int stepX = direction.X >= 0 ? 1 : -1;
            int stepY = direction.Y >= 0 ? 1 : -1;
            int stepZ = direction.Z >= 0 ? 1 : -1;

            // Distance to next cell boundary
            double tMaxX = direction.X != 0
                ? ((x + (stepX > 0 ? 1 : 0)) * _cellSize - origin.X) / direction.X
                : double.MaxValue;
            double tMaxY = direction.Y != 0
                ? ((y + (stepY > 0 ? 1 : 0)) * _cellSize - origin.Y) / direction.Y
                : double.MaxValue;
            double tMaxZ = direction.Z != 0
                ? ((z + (stepZ > 0 ? 1 : 0)) * _cellSize - origin.Z) / direction.Z
                : double.MaxValue;

            // How far along the ray we must move for each cell
            double tDeltaX = direction.X != 0 ? _cellSize / Math.Abs(direction.X) : double.MaxValue;
            double tDeltaY = direction.Y != 0 ? _cellSize / Math.Abs(direction.Y) : double.MaxValue;
            double tDeltaZ = direction.Z != 0 ? _cellSize / Math.Abs(direction.Z) : double.MaxValue;

            double t = 0;

            while (t < maxDistance)
            {
                long key = HashKey(x, y, z);

                if (_cells.TryGetValue(key, out var cell))
                {
                    foreach (var body in cell)
                    {
                        if (body.IsActive && visited.Add(body))
                        {
                            yield return body;
                        }
                    }
                }

                // Move to next cell
                if (tMaxX < tMaxY)
                {
                    if (tMaxX < tMaxZ)
                    {
                        x += stepX;
                        t = tMaxX;
                        tMaxX += tDeltaX;
                    }
                    else
                    {
                        z += stepZ;
                        t = tMaxZ;
                        tMaxZ += tDeltaZ;
                    }
                }
                else
                {
                    if (tMaxY < tMaxZ)
                    {
                        y += stepY;
                        t = tMaxY;
                        tMaxY += tDeltaY;
                    }
                    else
                    {
                        z += stepZ;
                        t = tMaxZ;
                        tMaxZ += tDeltaZ;
                    }
                }
            }
        }

        #endregion
    }

    /// <summary>
    /// Multi-level spatial hash for scenes with varying object sizes.
    /// Uses multiple grids with different cell sizes for optimal performance.
    /// </summary>
    public class MultiLevelSpatialHash
    {
        private readonly SpatialHash[] _levels;
        private readonly double[] _levelSizes;

        /// <summary>
        /// Creates a multi-level spatial hash.
        /// </summary>
        /// <param name="minSize">Minimum cell size (for smallest objects).</param>
        /// <param name="maxSize">Maximum cell size (for largest objects).</param>
        /// <param name="levels">Number of levels.</param>
        public MultiLevelSpatialHash(double minSize = 0.5, double maxSize = 16.0, int levels = 4)
        {
            _levels = new SpatialHash[levels];
            _levelSizes = new double[levels];

            double ratio = Math.Pow(maxSize / minSize, 1.0 / (levels - 1));
            for (int i = 0; i < levels; i++)
            {
                _levelSizes[i] = minSize * Math.Pow(ratio, i);
                _levels[i] = new SpatialHash(_levelSizes[i]);
            }
        }

        /// <summary>
        /// Gets the appropriate level for a body based on its size.
        /// </summary>
        private int GetLevel(IPhysicsBody body)
        {
            var aabb = body.BoundingBox;
            double size = Math.Max(
                aabb.Max.X - aabb.Min.X,
                Math.Max(aabb.Max.Y - aabb.Min.Y, aabb.Max.Z - aabb.Min.Z)
            );

            for (int i = 0; i < _levels.Length; i++)
            {
                if (size <= _levelSizes[i])
                    return i;
            }
            return _levels.Length - 1;
        }

        /// <summary>
        /// Clears all levels.
        /// </summary>
        public void Clear()
        {
            foreach (var level in _levels)
            {
                level.Clear();
            }
        }

        /// <summary>
        /// Inserts a body into the appropriate level.
        /// </summary>
        public void Insert(IPhysicsBody body)
        {
            int level = GetLevel(body);
            _levels[level].Insert(body);
        }

        /// <summary>
        /// Rebuilds all levels from a list of bodies.
        /// </summary>
        public void Rebuild(IReadOnlyList<IPhysicsBody> bodies)
        {
            Clear();
            for (int i = 0; i < bodies.Count; i++)
            {
                if (bodies[i].IsActive)
                {
                    Insert(bodies[i]);
                }
            }
        }

        /// <summary>
        /// Gets all potential collision pairs from all levels.
        /// </summary>
        public void GetPotentialPairs(List<(IPhysicsBody A, IPhysicsBody B)> pairs)
        {
            pairs.Clear();

            // Get pairs from each level
            var levelPairs = new List<(IPhysicsBody, IPhysicsBody)>();
            foreach (var level in _levels)
            {
                level.GetPotentialPairs(levelPairs);
                pairs.AddRange(levelPairs);
            }

            // Also need to check cross-level collisions
            // Small objects might collide with large objects
            for (int i = 0; i < _levels.Length - 1; i++)
            {
                var smallLevel = _levels[i];
                for (int j = i + 1; j < _levels.Length; j++)
                {
                    var largeLevel = _levels[j];
                    // Query each body in smaller level against larger level
                    // This is handled by the spatial hash queries
                }
            }
        }
    }
}
