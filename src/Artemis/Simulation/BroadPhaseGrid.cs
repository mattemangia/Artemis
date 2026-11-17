using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Simulation
{
    /// <summary>
    /// Uniform-grid broad phase shared by the 3D worlds.
    /// <para>
    /// Bodies are written into every cell their bounding box touches. Bodies that
    /// would cover a large number of cells (ground planes, walls) are kept in a
    /// separate list and tested against everything, which is far cheaper than
    /// flooding the grid. A pair that shares several cells is reported only from the
    /// first cell the two bodies have in common, so no hash set is needed to remove
    /// duplicates.
    /// </para>
    /// </summary>
    public sealed class BroadPhaseGrid
    {
        private readonly Dictionary<long, List<int>> _spatialGrid = new();
        private readonly Stack<List<int>> _cellListPool = new();
        private readonly List<int> _largeBodies = new();
        private readonly List<IPhysicsBody> _sweepCandidates = new();
        private readonly List<(int A, int B)> _pairs = new();
        private (int, int, int)[] _cellMin = new (int, int, int)[64];
        private bool[] _isLarge = new bool[64];
        private const double MaxCellsPerBody = 64;
        private const double CellCoordLimit = 1 << 20;
        private double _cellSize = 2.0;

        /// <summary>
        /// Edge length of the grid cells. Should be in the order of a typical body size.
        /// </summary>
        public double CellSize
        {
            get => _cellSize;
            set
            {
                if (!(value > 0) || double.IsInfinity(value))
                    throw new ArgumentOutOfRangeException(nameof(value), "Cell size must be positive and finite.");
                _cellSize = value;
            }
        }

        /// <summary>
        /// Candidate pairs (indices into the body list) found by the last <see cref="Build"/>.
        /// </summary>
        public IReadOnlyList<(int A, int B)> Pairs => _pairs;

        /// <summary>
        /// Number of occupied grid cells.
        /// </summary>
        public int OccupiedCellCount => _spatialGrid.Count;


        /// <summary>
        /// Bodies near a volume, taken from the grid of the last <see cref="Build"/>.
        /// The result may contain duplicates.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> CollectSweepCandidates(IReadOnlyList<IPhysicsBody> bodies, AABB swept)
        {
            double inv = 1.0 / CellSize;
            double minX = Math.Floor(swept.Min.X * inv), maxX = Math.Floor(swept.Max.X * inv);
            double minY = Math.Floor(swept.Min.Y * inv), maxY = Math.Floor(swept.Max.Y * inv);
            double minZ = Math.Floor(swept.Min.Z * inv), maxZ = Math.Floor(swept.Max.Z * inv);
            double cells = (maxX - minX + 1) * (maxY - minY + 1) * (maxZ - minZ + 1);

            // A sweep across a huge volume is cheaper to test against everything.
            if (!(cells <= 4096) || Math.Abs(minX) > CellCoordLimit || Math.Abs(maxX) > CellCoordLimit ||
                Math.Abs(minY) > CellCoordLimit || Math.Abs(maxY) > CellCoordLimit ||
                Math.Abs(minZ) > CellCoordLimit || Math.Abs(maxZ) > CellCoordLimit)
                return bodies;

            _sweepCandidates.Clear();
            foreach (int large in _largeBodies)
                _sweepCandidates.Add(bodies[large]);

            for (int cx = (int)minX; cx <= (int)maxX; cx++)
            {
                for (int cy = (int)minY; cy <= (int)maxY; cy++)
                {
                    for (int cz = (int)minZ; cz <= (int)maxZ; cz++)
                    {
                        if (!_spatialGrid.TryGetValue(GetCellKeyFromCoords(cx, cy, cz), out var list))
                            continue;

                        foreach (int index in list)
                            _sweepCandidates.Add(bodies[index]);
                    }
                }
            }

            return _sweepCandidates;
        }

        /// <summary>
        /// Rebuilds the grid and the list of candidate pairs. Each overlapping pair is
        /// reported once, with the lower body index first.
        /// </summary>
        public void Build(IReadOnlyList<IPhysicsBody> bodies)
        {
            foreach (var list in _spatialGrid.Values)
            {
                list.Clear();
                _cellListPool.Push(list);
            }
            _spatialGrid.Clear();
            _largeBodies.Clear();
            _pairs.Clear();

            int count = bodies.Count;
            if (_cellMin.Length < count)
            {
                _cellMin = new (int, int, int)[Math.Max(count, _cellMin.Length * 2)];
                _isLarge = new bool[_cellMin.Length];
            }
            Array.Clear(_isLarge, 0, count);

            double inv = 1.0 / CellSize;

            for (int i = 0; i < count; i++)
            {
                var body = bodies[i];
                if (!body.IsActive)
                    continue;

                var box = body.BoundingBox;
                double minX = Math.Floor(box.Min.X * inv), maxX = Math.Floor(box.Max.X * inv);
                double minY = Math.Floor(box.Min.Y * inv), maxY = Math.Floor(box.Max.Y * inv);
                double minZ = Math.Floor(box.Min.Z * inv), maxZ = Math.Floor(box.Max.Z * inv);

                double cells = (maxX - minX + 1) * (maxY - minY + 1) * (maxZ - minZ + 1);

                // Bodies that span many cells (ground planes, walls) or whose bounds are
                // not finite are tested against everything instead of flooding the grid.
                if (!(cells <= MaxCellsPerBody) || Math.Abs(minX) > CellCoordLimit || Math.Abs(maxX) > CellCoordLimit ||
                    Math.Abs(minY) > CellCoordLimit || Math.Abs(maxY) > CellCoordLimit ||
                    Math.Abs(minZ) > CellCoordLimit || Math.Abs(maxZ) > CellCoordLimit)
                {
                    if (double.IsNaN(cells))
                        continue;
                    _largeBodies.Add(i);
                    _isLarge[i] = true;
                    continue;
                }

                var min = ((int)minX, (int)minY, (int)minZ);
                var max = ((int)maxX, (int)maxY, (int)maxZ);
                _cellMin[i] = min;

                for (int cx = min.Item1; cx <= max.Item1; cx++)
                {
                    for (int cy = min.Item2; cy <= max.Item2; cy++)
                    {
                        for (int cz = min.Item3; cz <= max.Item3; cz++)
                        {
                            long key = GetCellKeyFromCoords(cx, cy, cz);
                            if (!_spatialGrid.TryGetValue(key, out var list))
                            {
                                list = _cellListPool.Count > 0 ? _cellListPool.Pop() : new List<int>(8);
                                _spatialGrid[key] = list;
                            }
                            list.Add(i);
                        }
                    }
                }
            }

            // Pairs sharing a cell. Two bodies can share several cells, so a pair is
            // reported only from the first cell their cell ranges have in common.
            foreach (var cell in _spatialGrid)
            {
                var indices = cell.Value;
                if (indices.Count < 2)
                    continue;

                long cellKey = cell.Key;

                for (int a = 0; a < indices.Count; a++)
                {
                    int i = indices[a];
                    for (int b = a + 1; b < indices.Count; b++)
                    {
                        int j = indices[b];

                        var minI = _cellMin[i];
                        var minJ = _cellMin[j];
                        long firstShared = GetCellKeyFromCoords(
                            Math.Max(minI.Item1, minJ.Item1),
                            Math.Max(minI.Item2, minJ.Item2),
                            Math.Max(minI.Item3, minJ.Item3));
                        if (firstShared != cellKey)
                            continue;

                        AddPairIfOverlapping(bodies, i, j);
                    }
                }
            }

            // Large bodies against everything else.
            for (int a = 0; a < _largeBodies.Count; a++)
            {
                int large = _largeBodies[a];
                for (int k = 0; k < count; k++)
                {
                    if (k == large)
                        continue;

                    // Large-large pairs are visited from both sides: keep one.
                    if (k < large && _isLarge[k])
                        continue;

                    AddPairIfOverlapping(bodies, large, k);
                }
            }
        }

        private void AddPairIfOverlapping(IReadOnlyList<IPhysicsBody> bodies, int i, int j)
        {
            var bodyA = bodies[i];
            var bodyB = bodies[j];

            if (!bodyA.IsActive || !bodyB.IsActive)
                return;

            // At least one side must be able to react.
            bool movingA = bodyA.BodyType == BodyType.Dynamic && !bodyA.IsSleeping;
            bool movingB = bodyB.BodyType == BodyType.Dynamic && !bodyB.IsSleeping;
            if (!movingA && !movingB)
            {
                // A moving kinematic body must still be able to wake a sleeper.
                bool sleeperAndMover =
                    (bodyA.BodyType == BodyType.Dynamic && SleepManager.IsMovingKinematic(bodyB)) ||
                    (bodyB.BodyType == BodyType.Dynamic && SleepManager.IsMovingKinematic(bodyA));
                if (!sleeperAndMover)
                    return;
            }

            if (!bodyA.BoundingBox.Intersects(bodyB.BoundingBox))
                return;

            _pairs.Add(i < j ? (i, j) : (j, i));
        }

        private static long GetCellKeyFromCoords(int cx, int cy, int cz)
        {
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (long)(cz & 0x1FFFFF);
        }
    }
}
