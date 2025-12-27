using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Artemis.Bodies;
using Artemis.Collision;

namespace Artemis.Simulation
{
    /// <summary>
    /// Island solver that groups connected bodies for independent parallel solving.
    /// Bodies in the same island share constraints (contacts, joints) and must be solved together.
    /// Different islands can be solved in parallel.
    /// </summary>
    public class IslandSolver
    {
        #region Nested Types

        /// <summary>
        /// Represents a group of connected bodies that must be solved together.
        /// </summary>
        public class Island
        {
            /// <summary>Bodies in this island.</summary>
            public readonly List<IPhysicsBody> Bodies = new();

            /// <summary>Collision pairs in this island.</summary>
            public readonly List<CollisionInfo> Collisions = new();

            /// <summary>Whether all bodies are sleeping.</summary>
            public bool IsSleeping { get; internal set; }

            /// <summary>
            /// Clears the island for reuse.
            /// </summary>
            public void Clear()
            {
                Bodies.Clear();
                Collisions.Clear();
                IsSleeping = false;
            }
        }

        #endregion

        #region Fields

        private readonly List<Island> _islands;
        private readonly List<Island> _islandPool;
        private readonly Dictionary<IPhysicsBody, int> _bodyToIsland;
        private readonly int[] _parent;
        private readonly int[] _rank;
        private int _bodyCount;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the list of islands.
        /// </summary>
        public IReadOnlyList<Island> Islands => _islands;

        /// <summary>
        /// Gets the number of active islands.
        /// </summary>
        public int IslandCount => _islands.Count;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new island solver with a maximum body capacity.
        /// </summary>
        /// <param name="maxBodies">Maximum number of bodies.</param>
        public IslandSolver(int maxBodies = 10000)
        {
            _islands = new List<Island>(64);
            _islandPool = new List<Island>(64);
            _bodyToIsland = new Dictionary<IPhysicsBody, int>(maxBodies);
            _parent = new int[maxBodies];
            _rank = new int[maxBodies];
        }

        #endregion

        #region Union-Find (Disjoint Set)

        /// <summary>
        /// Finds the root of the set containing element i with path compression.
        /// </summary>
        private int Find(int i)
        {
            if (_parent[i] != i)
            {
                _parent[i] = Find(_parent[i]); // Path compression
            }
            return _parent[i];
        }

        /// <summary>
        /// Unites the sets containing elements i and j using union by rank.
        /// </summary>
        private void Union(int i, int j)
        {
            int rootI = Find(i);
            int rootJ = Find(j);

            if (rootI != rootJ)
            {
                // Union by rank
                if (_rank[rootI] < _rank[rootJ])
                {
                    _parent[rootI] = rootJ;
                }
                else if (_rank[rootI] > _rank[rootJ])
                {
                    _parent[rootJ] = rootI;
                }
                else
                {
                    _parent[rootJ] = rootI;
                    _rank[rootI]++;
                }
            }
        }

        #endregion

        #region Island Building

        /// <summary>
        /// Builds islands from bodies and their collision pairs.
        /// </summary>
        /// <param name="bodies">All bodies in the simulation.</param>
        /// <param name="collisions">All collision pairs.</param>
        public void BuildIslands(
            IReadOnlyList<IPhysicsBody> bodies,
            IReadOnlyList<CollisionInfo> collisions)
        {
            // Return islands to pool
            foreach (var island in _islands)
            {
                island.Clear();
                _islandPool.Add(island);
            }
            _islands.Clear();
            _bodyToIsland.Clear();

            _bodyCount = bodies.Count;
            if (_bodyCount == 0) return;

            // Initialize union-find
            for (int i = 0; i < _bodyCount; i++)
            {
                _parent[i] = i;
                _rank[i] = 0;
                _bodyToIsland[bodies[i]] = i;
            }

            // Union connected bodies through collisions
            foreach (var collision in collisions)
            {
                if (_bodyToIsland.TryGetValue(collision.BodyA, out int indexA) &&
                    _bodyToIsland.TryGetValue(collision.BodyB, out int indexB))
                {
                    // Don't union static bodies - they belong to their own island
                    bool aStatic = collision.BodyA.BodyType == BodyType.Static;
                    bool bStatic = collision.BodyB.BodyType == BodyType.Static;

                    if (!aStatic && !bStatic)
                    {
                        Union(indexA, indexB);
                    }
                    else if (!aStatic)
                    {
                        // A is dynamic, B is static - A stays in its island
                    }
                    else if (!bStatic)
                    {
                        // B is dynamic, A is static - B stays in its island
                    }
                }
            }

            // Build islands from union-find results
            var rootToIsland = new Dictionary<int, Island>(_bodyCount);

            for (int i = 0; i < _bodyCount; i++)
            {
                var body = bodies[i];
                if (!body.IsActive) continue;
                if (body.BodyType == BodyType.Static) continue; // Skip static bodies for islands

                int root = Find(i);

                if (!rootToIsland.TryGetValue(root, out var island))
                {
                    island = GetOrCreateIsland();
                    rootToIsland[root] = island;
                    _islands.Add(island);
                }

                island.Bodies.Add(body);
            }

            // Assign collisions to islands
            foreach (var collision in collisions)
            {
                // Determine which island this collision belongs to
                int islandIndex = -1;

                if (collision.BodyA.BodyType != BodyType.Static &&
                    _bodyToIsland.TryGetValue(collision.BodyA, out int indexA))
                {
                    int root = Find(indexA);
                    if (rootToIsland.TryGetValue(root, out var island))
                    {
                        islandIndex = _islands.IndexOf(island);
                        island.Collisions.Add(collision);
                    }
                }
                else if (collision.BodyB.BodyType != BodyType.Static &&
                         _bodyToIsland.TryGetValue(collision.BodyB, out int indexB))
                {
                    int root = Find(indexB);
                    if (rootToIsland.TryGetValue(root, out var island))
                    {
                        if (islandIndex == -1) // Not already added
                        {
                            island.Collisions.Add(collision);
                        }
                    }
                }
            }

            // Determine sleeping status for each island
            foreach (var island in _islands)
            {
                island.IsSleeping = true;
                foreach (var body in island.Bodies)
                {
                    if (!body.IsSleeping)
                    {
                        island.IsSleeping = false;
                        break;
                    }
                }
            }
        }

        private Island GetOrCreateIsland()
        {
            if (_islandPool.Count > 0)
            {
                var island = _islandPool[_islandPool.Count - 1];
                _islandPool.RemoveAt(_islandPool.Count - 1);
                return island;
            }
            return new Island();
        }

        #endregion

        #region Parallel Solving

        /// <summary>
        /// Solves all islands in parallel.
        /// </summary>
        /// <param name="solveAction">Action to execute for each island.</param>
        public void SolveParallel(Action<Island> solveAction)
        {
            // Filter out sleeping islands
            var activeIslands = new List<Island>();
            foreach (var island in _islands)
            {
                if (!island.IsSleeping)
                {
                    activeIslands.Add(island);
                }
            }

            if (activeIslands.Count == 0) return;

            // Solve islands in parallel
            Parallel.ForEach(activeIslands, island =>
            {
                solveAction(island);
            });
        }

        /// <summary>
        /// Solves all islands in parallel with custom parallelism.
        /// </summary>
        public void SolveParallel(Action<Island> solveAction, int maxDegreeOfParallelism)
        {
            var activeIslands = new List<Island>();
            foreach (var island in _islands)
            {
                if (!island.IsSleeping)
                {
                    activeIslands.Add(island);
                }
            }

            if (activeIslands.Count == 0) return;

            var options = new ParallelOptions
            {
                MaxDegreeOfParallelism = maxDegreeOfParallelism
            };

            Parallel.ForEach(activeIslands, options, island =>
            {
                solveAction(island);
            });
        }

        /// <summary>
        /// Solves all islands sequentially (for debugging or single-threaded mode).
        /// </summary>
        public void SolveSequential(Action<Island> solveAction)
        {
            foreach (var island in _islands)
            {
                if (!island.IsSleeping)
                {
                    solveAction(island);
                }
            }
        }

        #endregion

        #region Statistics

        /// <summary>
        /// Gets the average island size.
        /// </summary>
        public double AverageIslandSize
        {
            get
            {
                if (_islands.Count == 0) return 0;
                int total = 0;
                foreach (var island in _islands)
                {
                    total += island.Bodies.Count;
                }
                return (double)total / _islands.Count;
            }
        }

        /// <summary>
        /// Gets the largest island size.
        /// </summary>
        public int LargestIslandSize
        {
            get
            {
                int max = 0;
                foreach (var island in _islands)
                {
                    if (island.Bodies.Count > max)
                        max = island.Bodies.Count;
                }
                return max;
            }
        }

        /// <summary>
        /// Gets the number of sleeping islands.
        /// </summary>
        public int SleepingIslandCount
        {
            get
            {
                int count = 0;
                foreach (var island in _islands)
                {
                    if (island.IsSleeping) count++;
                }
                return count;
            }
        }

        #endregion
    }
}
