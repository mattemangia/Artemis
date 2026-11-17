using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Compute;

namespace Artemis.Physics2D
{
    /// <summary>
    /// High-performance 2D physics world with multi-threading, SIMD, and GPU support.
    /// </summary>
    public class PhysicsWorld2D
    {
        #region Fields

        private readonly List<RigidBody2D> _bodies = new();
        private readonly List<Joint2D> _joints = new();
        private readonly List<AreaEffector2D> _areaEffectors = new();
        private readonly List<Manifold2D> _manifolds = new();
        private readonly Dictionary<(string, string), Manifold2D> _manifoldCache = new();

        // Spatial hash for broad phase
        private readonly Dictionary<long, List<int>> _grid = new();
        private readonly Stack<List<int>> _cellListPool = new();
        private readonly List<int> _largeBodies = new();
        private readonly List<(int A, int B)> _pairs = new();
        private (int, int)[] _cellMin = new (int, int)[64];
        private Manifold2D?[] _narrowResults = new Manifold2D?[64];
        private const double MaxCellsPerBody = 1024;
        private const double CellCoordLimit = 1 << 30;
        private const int ParallelThreshold = 256;
        private const double ContactMatchDistance = 0.05;

        // Island based sleeping
        private int[] _islandParent = new int[64];
        private byte[] _islandFlag = new byte[64];
        private float _cellSize = 2.0f;
        private float _invCellSize = 0.5f;

        // Thread synchronization
        private readonly ReaderWriterLockSlim _bodiesLock = new();

        // GPU compute (optional)

        #endregion

        #region Properties

        /// <summary>Gravity vector applied to all dynamic bodies.</summary>
        public Vector2D Gravity { get; set; } = new(0, -9.81);

        /// <summary>Number of velocity iterations per step.</summary>
        public int VelocityIterations { get; set; } = 8;

        /// <summary>Number of position iterations per step.</summary>
        public int PositionIterations { get; set; } = 3;

        /// <summary>All bodies in the world (thread-safe copy).</summary>
        public IReadOnlyList<RigidBody2D> Bodies
        {
            get
            {
                if (_bodiesLock.IsReadLockHeld || _bodiesLock.IsWriteLockHeld)
                    return _bodies.ToArray();

                _bodiesLock.EnterReadLock();
                try { return _bodies.ToArray(); }
                finally { _bodiesLock.ExitReadLock(); }
            }
        }

        /// <summary>All joints in the world.</summary>
        public IReadOnlyList<Joint2D> Joints => _joints;

        /// <summary>All area effectors in the world.</summary>
        public IReadOnlyList<AreaEffector2D> AreaEffectors => _areaEffectors;

        /// <summary>All active collision manifolds.</summary>
        public IReadOnlyList<Manifold2D> Manifolds => _manifolds;

        /// <summary>Whether to allow bodies to sleep.</summary>
        public bool AllowSleep { get; set; } = true;

        /// <summary>Total simulation time elapsed.</summary>
        public double SimulationTime { get; private set; }

        /// <summary>
        /// Kept for compatibility. The 2D world has no GPU path and always runs on the
        /// CPU, so the value stays false.
        /// </summary>
        public bool UseGPU
        {
            get => false;
            set { }
        }

        /// <summary>Whether to use continuous collision detection.</summary>
        public bool UseCCD { get; set; } = true;

        /// <summary>
        /// Approach speed (m/s) below which a contact is treated as inelastic, so
        /// resting bodies do not jitter.
        /// </summary>
        public double RestitutionThreshold { get; set; } = 1.0;

        /// <summary>
        /// Overlap (m) tolerated without positional correction.
        /// </summary>
        public double PenetrationSlop { get; set; } = 0.005;

        /// <summary>
        /// Linear speed (m/s) below which a body counts as resting.
        /// </summary>
        public double SleepVelocityThreshold { get; set; } = 0.05;

        /// <summary>
        /// Angular speed (rad/s) below which a body counts as resting.
        /// </summary>
        public double SleepAngularVelocityThreshold { get; set; } = 0.1;

        /// <summary>
        /// Time (s) a whole island has to rest before it is put to sleep.
        /// </summary>
        public double SleepTimeThreshold { get; set; } = 0.5;

        private bool _useAdvancedSolver;

        /// <summary>Enable advanced solver settings for higher stability.</summary>
        public bool UseAdvancedSolver
        {
            get => _useAdvancedSolver;
            set
            {
                _useAdvancedSolver = value;
                if (value)
                {
                    VelocityIterations = 12;
                    PositionIterations = 6;
                }
                else
                {
                    VelocityIterations = 8;
                    PositionIterations = 3;
                }
            }
        }

        #endregion

        #region Events

        /// <summary>Raised when a collision begins.</summary>
        public event EventHandler<CollisionEvent2D>? CollisionBegin;

        /// <summary>Raised when a collision ends.</summary>
        public event EventHandler<CollisionEvent2D>? CollisionEnd;

        /// <summary>Raised during collision (for modifying contact data).</summary>
        public event EventHandler<CollisionEvent2D>? PreSolve;

        /// <summary>Raised after collision is resolved.</summary>
        public event EventHandler<CollisionEvent2D>? PostSolve;

        /// <summary>Raised when a trigger is entered.</summary>
        public event EventHandler<CollisionEvent2D>? TriggerEnter;

        /// <summary>Raised when a trigger is exited.</summary>
        public event EventHandler<CollisionEvent2D>? TriggerExit;

        /// <summary>Compatibility event for collision begin.</summary>
        public event EventHandler<CollisionEvent2D>? OnCollisionEnter
        {
            add => CollisionBegin += value;
            remove => CollisionBegin -= value;
        }

        /// <summary>Compatibility event for collision end.</summary>
        public event EventHandler<CollisionEvent2D>? OnCollisionExit
        {
            add => CollisionEnd += value;
            remove => CollisionEnd -= value;
        }

        /// <summary>Compatibility event for pre-solve.</summary>
        public event EventHandler<CollisionEvent2D>? OnPreSolve
        {
            add => PreSolve += value;
            remove => PreSolve -= value;
        }

        /// <summary>Compatibility event for post-solve.</summary>
        public event EventHandler<CollisionEvent2D>? OnPostSolve
        {
            add => PostSolve += value;
            remove => PostSolve -= value;
        }

        /// <summary>Compatibility event for trigger enter.</summary>
        public event EventHandler<CollisionEvent2D>? OnTriggerEnter
        {
            add => TriggerEnter += value;
            remove => TriggerEnter -= value;
        }

        /// <summary>Compatibility event for trigger exit.</summary>
        public event EventHandler<CollisionEvent2D>? OnTriggerExit
        {
            add => TriggerExit += value;
            remove => TriggerExit -= value;
        }

        #endregion

        #region Constructor

        public PhysicsWorld2D()
        {
        }

        public PhysicsWorld2D(Vector2D gravity) : this()
        {
            Gravity = gravity;
        }

        #endregion

        #region Body Management

        /// <summary>
        /// Adds a body to the world.
        /// </summary>
        public void AddBody(RigidBody2D body)
        {
            _bodiesLock.EnterWriteLock();
            try
            {
                if (!_bodies.Contains(body))
                    _bodies.Add(body);
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        /// <summary>
        /// Removes a body from the world.
        /// </summary>
        public bool RemoveBody(RigidBody2D body)
        {
            _bodiesLock.EnterWriteLock();
            try
            {
                _joints.RemoveAll(j => j.BodyA == body || j.BodyB == body);
                foreach (var manifold in _manifolds)
                {
                    if (manifold.BodyA == body || manifold.BodyB == body)
                        _manifoldCache.Remove(GetManifoldKey(manifold.BodyA, manifold.BodyB));
                }
                _manifolds.RemoveAll(m => m.BodyA == body || m.BodyB == body);

                bool removed = _bodies.Remove(body);
                if (removed)
                {
                    // Whatever was resting on or hanging from the body lost its support.
                    foreach (var other in _bodies)
                    {
                        if (other.IsSleeping)
                            other.WakeUp();
                    }
                }
                return removed;
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        /// <summary>
        /// Gets a body by ID.
        /// </summary>
        public RigidBody2D? GetBody(string id)
        {
            _bodiesLock.EnterReadLock();
            try { return _bodies.Find(b => b.Id == id); }
            finally { _bodiesLock.ExitReadLock(); }
        }

        /// <summary>
        /// Clears all bodies from the world.
        /// </summary>
        public void ClearBodies()
        {
            _bodiesLock.EnterWriteLock();
            try
            {
                _bodies.Clear();
                _joints.Clear();
                _manifolds.Clear();
                _manifoldCache.Clear();
            }
            finally { _bodiesLock.ExitWriteLock(); }
        }

        #endregion

        #region Joint Management

        public void AddJoint(Joint2D joint)
        {
            if (!_joints.Contains(joint))
                _joints.Add(joint);
        }

        public bool RemoveJoint(Joint2D joint) => _joints.Remove(joint);

        public DistanceJoint2D CreateDistanceJoint(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchorA, Vector2D anchorB)
        {
            var joint = new DistanceJoint2D(bodyA, bodyB, anchorA, anchorB);
            AddJoint(joint);
            return joint;
        }

        public RevoluteJoint2D CreateRevoluteJoint(RigidBody2D bodyA, RigidBody2D bodyB, Vector2D anchor)
        {
            var joint = new RevoluteJoint2D(bodyA, bodyB, anchor);
            AddJoint(joint);
            return joint;
        }

        #endregion

        #region Area Effector Management

        public void AddAreaEffector(AreaEffector2D effector)
        {
            if (!_areaEffectors.Contains(effector))
                _areaEffectors.Add(effector);
        }

        public bool RemoveAreaEffector(AreaEffector2D effector) => _areaEffectors.Remove(effector);

        #endregion

        #region Simulation

        /// <summary>
        /// Advances the simulation by one step.
        /// <para>
        /// Pipeline: forces -> velocities -> collision detection -> joint and contact
        /// solver -> positions (with continuous collision detection for bullets) ->
        /// penetration recovery -> sleeping.
        /// </para>
        /// </summary>
        public void Step(double dt)
        {
            if (dt <= 0) return;

            _bodiesLock.EnterReadLock();
            try
            {
                for (int i = 0; i < _bodies.Count; i++)
                    _bodies[i].WorldIndex = i;

                // Phase 1: forces and velocities
                ApplyAreaEffectorsParallel(dt);
                IntegrateVelocitiesParallel(dt);

                // Phase 2: collision detection
                BroadPhase();
                NarrowPhase();

                // Phase 3: wake whatever is being disturbed
                BuildIslands();
                WakeDisturbedIslands();

                // Phase 4: prepare constraints
                // Contacts are prepared before anything is warm started: restitution
                // depends on the approach speed, which the warm-start impulses of
                // neighbouring contacts (and joints) would otherwise distort.
                foreach (var manifold in _manifolds)
                {
                    if (manifold.IsActive && !manifold.IsTrigger)
                        ContactSolver2D.Prepare(manifold, RestitutionThreshold);
                }

                foreach (var joint in _joints)
                {
                    if (joint.IsActive)
                        joint.PreSolve(dt);
                }

                foreach (var manifold in _manifolds)
                {
                    if (manifold.IsActive && !manifold.IsTrigger)
                        ContactSolver2D.WarmStart(manifold);
                }

                // Phase 5: velocity iterations
                for (int i = 0; i < VelocityIterations; i++)
                {
                    foreach (var joint in _joints)
                    {
                        if (joint.IsActive)
                            joint.SolveVelocity();
                    }

                    foreach (var manifold in _manifolds)
                    {
                        if (manifold.IsActive && !manifold.IsTrigger)
                            ContactSolver2D.SolveVelocity(manifold);
                    }
                }

                // Phase 6: positions
                IntegratePositions(dt);

                // Phase 7: penetration recovery and joint drift
                CorrectContactPositions();

                for (int i = 0; i < PositionIterations; i++)
                {
                    bool jointsSolved = true;
                    foreach (var joint in _joints)
                    {
                        if (joint.IsActive)
                            jointsSolved &= joint.SolvePosition();
                    }

                    if (jointsSolved) break;
                }

                // Phase 8: Post-solve events
                if (PostSolve != null)
                {
                    foreach (var manifold in _manifolds)
                    {
                        if (manifold.IsActive && !manifold.IsTrigger)
                            PostSolve.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                    }
                }

                // Phase 9: sleeping
                UpdateSleeping(dt);
            }
            finally { _bodiesLock.ExitReadLock(); }

            SimulationTime += dt;
        }

        /// <summary>
        /// Steps the simulation with fixed timestep and accumulator.
        /// </summary>
        public void StepFixed(double dt, double fixedDt = 1.0 / 60.0)
        {
            while (dt >= fixedDt)
            {
                Step(fixedDt);
                dt -= fixedDt;
            }

            if (dt > PhysicsConstants.Epsilon)
                Step(dt);
        }

        #endregion

        #region Broad Phase

        private void BroadPhase()
        {
            foreach (var list in _grid.Values)
            {
                list.Clear();
                _cellListPool.Push(list);
            }
            _grid.Clear();
            _largeBodies.Clear();
            _pairs.Clear();

            int count = _bodies.Count;
            if (_cellMin.Length < count)
                _cellMin = new (int, int)[Math.Max(count, _cellMin.Length * 2)];

            double inv = _invCellSize;

            for (int i = 0; i < count; i++)
            {
                var body = _bodies[i];
                if (!body.IsActive || body.Shape == null) continue;

                var box = body.AABB;
                double minX = Math.Floor(box.Min.X * inv), maxX = Math.Floor(box.Max.X * inv);
                double minY = Math.Floor(box.Min.Y * inv), maxY = Math.Floor(box.Max.Y * inv);
                double cells = (maxX - minX + 1) * (maxY - minY + 1);

                if (double.IsNaN(cells))
                    continue;

                // Very large bodies (long floors, chains) are tested against everything
                // instead of being written into thousands of cells.
                if (cells > MaxCellsPerBody || Math.Abs(minX) > CellCoordLimit || Math.Abs(maxX) > CellCoordLimit ||
                    Math.Abs(minY) > CellCoordLimit || Math.Abs(maxY) > CellCoordLimit)
                {
                    _largeBodies.Add(i);
                    continue;
                }

                _cellMin[i] = ((int)minX, (int)minY);

                for (int x = (int)minX; x <= (int)maxX; x++)
                {
                    for (int y = (int)minY; y <= (int)maxY; y++)
                    {
                        long key = GetCellKey(x, y);
                        if (!_grid.TryGetValue(key, out var list))
                        {
                            list = _cellListPool.Count > 0 ? _cellListPool.Pop() : new List<int>(8);
                            _grid[key] = list;
                        }
                        list.Add(i);
                    }
                }
            }

            // Two bodies can share several cells: the pair is taken only from the first
            // cell their ranges have in common, so each pair is produced exactly once.
            foreach (var cell in _grid)
            {
                var indices = cell.Value;
                for (int a = 0; a < indices.Count; a++)
                {
                    int i = indices[a];
                    for (int b = a + 1; b < indices.Count; b++)
                    {
                        int j = indices[b];
                        long firstShared = GetCellKey(
                            Math.Max(_cellMin[i].Item1, _cellMin[j].Item1),
                            Math.Max(_cellMin[i].Item2, _cellMin[j].Item2));
                        if (firstShared == cell.Key)
                            AddPairIfOverlapping(i, j);
                    }
                }
            }

            for (int a = 0; a < _largeBodies.Count; a++)
            {
                int large = _largeBodies[a];
                for (int k = 0; k < count; k++)
                {
                    if (k == large) continue;
                    if (k < large && _largeBodies.Contains(k)) continue;
                    AddPairIfOverlapping(large, k);
                }
            }

            // Fixed order regardless of how the grid enumerates.
            _pairs.Sort(static (x, y) => x.A != y.A ? x.A.CompareTo(y.A) : x.B.CompareTo(y.B));
        }

        private void AddPairIfOverlapping(int i, int j)
        {
            var a = _bodies[i];
            var b = _bodies[j];

            if (!a.IsActive || !b.IsActive || a.Shape == null || b.Shape == null) return;
            if (a.BodyType != BodyType2D.Dynamic && b.BodyType != BodyType2D.Dynamic)
            {
                // Static/kinematic pairs only matter for trigger events.
                if (!a.IsTrigger && !b.IsTrigger) return;
            }
            if (!a.ShouldCollide(b)) return;
            if (!a.AABB.Overlaps(b.AABB)) return;
            if (!ShouldJointedBodiesCollide(a, b)) return;

            _pairs.Add(i < j ? (i, j) : (j, i));
        }

        private bool ShouldJointedBodiesCollide(RigidBody2D a, RigidBody2D b)
        {
            foreach (var joint in _joints)
            {
                if (joint.CollideConnected || !joint.IsActive) continue;
                if ((joint.BodyA == a && joint.BodyB == b) || (joint.BodyA == b && joint.BodyB == a))
                    return false;
            }
            return true;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static long GetCellKey(int cx, int cy)
        {
            unchecked
            {
                return ((long)(uint)cx << 32) | (uint)cy;
            }
        }

        #endregion

        #region Narrow Phase

        private void NarrowPhase()
        {
            foreach (var manifold in _manifolds)
                manifold.IsActive = false;

            int pairCount = _pairs.Count;
            if (_narrowResults.Length < pairCount)
                _narrowResults = new Manifold2D?[Math.Max(pairCount, _narrowResults.Length * 2)];

            var results = _narrowResults;

            void Detect(int p)
            {
                var (i, j) = _pairs[p];
                var manifold = new Manifold2D();
                results[p] = CollisionDetector2D.DetectCollision(_bodies[i], _bodies[j], manifold) ? manifold : null;
            }

            // Results are stored by pair index, so the outcome does not depend on
            // thread scheduling.
            if (pairCount >= ParallelThreshold)
                Parallel.For(0, pairCount, Detect);
            else
                for (int p = 0; p < pairCount; p++) Detect(p);

            for (int p = 0; p < pairCount; p++)
            {
                var manifold = results[p];
                results[p] = null;
                if (manifold == null)
                    continue;

                manifold.IsActive = true;
                manifold.IsTrigger = manifold.BodyA.IsTrigger || manifold.BodyB.IsTrigger;

                var a = manifold.BodyA;
                var b = manifold.BodyB;
                for (int i = 0; i < manifold.ContactCount; i++)
                {
                    manifold.Contacts[i].LocalA = a.WorldToLocal(manifold.Contacts[i].Point);
                    manifold.Contacts[i].LocalB = b.WorldToLocal(manifold.Contacts[i].Point);
                }

                var key = GetManifoldKey(a, b);

                if (!_manifoldCache.TryGetValue(key, out var cached))
                {
                    _manifoldCache[key] = manifold;
                    _manifolds.Add(manifold);
                    cached = manifold;
                    if (manifold.IsTrigger)
                        TriggerEnter?.Invoke(this, new CollisionEvent2D(a, b, manifold));
                    else
                        CollisionBegin?.Invoke(this, new CollisionEvent2D(a, b, manifold));
                }
                else
                {
                    // Carry the accumulated impulses over to the matching new points
                    // (warm starting). Points are matched in the bodies' local frames.
                    for (int i = 0; i < manifold.ContactCount; i++)
                    {
                        ref var fresh = ref manifold.Contacts[i];
                        double best = ContactMatchDistance * ContactMatchDistance;
                        for (int k = 0; k < cached.ContactCount; k++)
                        {
                            ref var old = ref cached.Contacts[k];
                            double d = Math.Min((old.LocalA - fresh.LocalA).MagnitudeSquared,
                                                (old.LocalB - fresh.LocalB).MagnitudeSquared);
                            if (d < best)
                            {
                                best = d;
                                fresh.NormalImpulse = old.NormalImpulse;
                                fresh.TangentImpulse = old.TangentImpulse;
                            }
                        }
                    }

                    cached.BodyA = a;
                    cached.BodyB = b;
                    cached.Normal = manifold.Normal;
                    cached.Penetration = manifold.Penetration;
                    cached.ContactCount = manifold.ContactCount;
                    for (int i = 0; i < manifold.ContactCount; i++)
                        cached.Contacts[i] = manifold.Contacts[i];
                    cached.Friction = manifold.Friction;
                    cached.Restitution = manifold.Restitution;
                    cached.IsActive = true;
                    cached.IsTrigger = manifold.IsTrigger;
                }

                if (!cached.IsTrigger)
                    PreSolve?.Invoke(this, new CollisionEvent2D(a, b, cached));
            }

            // Clean up ended collisions
            for (int i = _manifolds.Count - 1; i >= 0; i--)
            {
                var manifold = _manifolds[i];
                if (!manifold.IsActive)
                {
                    if (manifold.IsTrigger)
                        TriggerExit?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                    else
                        CollisionEnd?.Invoke(this, new CollisionEvent2D(manifold.BodyA, manifold.BodyB, manifold));
                    _manifoldCache.Remove(GetManifoldKey(manifold.BodyA, manifold.BodyB));
                    _manifolds.RemoveAt(i);
                }
            }
        }

        private static (string, string) GetManifoldKey(RigidBody2D a, RigidBody2D b)
        {
            return string.CompareOrdinal(a.Id, b.Id) < 0 ? (a.Id, b.Id) : (b.Id, a.Id);
        }

        #endregion

        #region Integration

        private void ApplyAreaEffectorsParallel(double dt)
        {
            if (_areaEffectors.Count == 0) return;

            var enabledEffectors = _areaEffectors.Where(e => e.Enabled).ToArray();
            var dynamicBodies = _bodies.Where(b => b.BodyType == BodyType2D.Dynamic && b.IsActive && !b.IsSleeping).ToArray();

            Parallel.ForEach(dynamicBodies, body =>
            {
                foreach (var effector in enabledEffectors)
                {
                    if (effector.AffectsBody(body))
                    {
                        effector.ApplyForce(body, dt);
                    }
                }
            });
        }

        private void IntegrateVelocitiesParallel(double dt)
        {
            var gravity = Gravity;

            if (_bodies.Count >= ParallelThreshold)
                Parallel.ForEach(_bodies, body => body.IntegrateVelocity(dt, gravity));
            else
                foreach (var body in _bodies) body.IntegrateVelocity(dt, gravity);
        }

        private void IntegratePositions(double dt)
        {
            bool ccd = UseCCD;

            foreach (var body in _bodies)
            {
                if (!body.IsActive) continue;

                // Bullets are swept against everything. Any other fast body is still
                // swept against static and kinematic geometry, so debris from an impact
                // cannot be knocked straight through the floor.
                double stepTime = dt;
                if (ccd && body.BodyType == BodyType2D.Dynamic && !body.IsSleeping)
                    stepTime = ContinuousCollisionDetector2D.TimeOfImpact(body, _bodies, dt,
                        againstDynamicBodies: body.IsBullet || body.UseCCD);

                body.IntegratePosition(stepTime);
            }
        }

        private void CorrectContactPositions()
        {
            if (_manifolds.Count == 0) return;

            foreach (var manifold in _manifolds)
            {
                manifold.BodyA.CorrectionLinear = Vector2D.Zero;
                manifold.BodyA.CorrectionAngular = 0;
                manifold.BodyB.CorrectionLinear = Vector2D.Zero;
                manifold.BodyB.CorrectionAngular = 0;
            }

            for (int i = 0; i < PositionIterations; i++)
            {
                foreach (var manifold in _manifolds)
                {
                    if (manifold.IsActive && !manifold.IsTrigger)
                        ContactSolver2D.SolvePosition(manifold, PenetrationSlop, 0.3, 0.25);
                }
            }

            foreach (var manifold in _manifolds)
            {
                ApplyCorrection(manifold.BodyA);
                ApplyCorrection(manifold.BodyB);
            }
        }

        private static void ApplyCorrection(RigidBody2D body)
        {
            if (body.CorrectionLinear.MagnitudeSquared <= 0 && body.CorrectionAngular == 0)
                return;

            body.ApplyPositionCorrection(body.CorrectionLinear, body.CorrectionAngular);
            body.CorrectionLinear = Vector2D.Zero;
            body.CorrectionAngular = 0;
        }

        #endregion

        #region Islands and Sleeping

        private int FindRoot(int index)
        {
            while (_islandParent[index] != index)
            {
                _islandParent[index] = _islandParent[_islandParent[index]];
                index = _islandParent[index];
            }
            return index;
        }

        private void Union(RigidBody2D a, RigidBody2D? b)
        {
            if (b == null || a.BodyType != BodyType2D.Dynamic || b.BodyType != BodyType2D.Dynamic)
                return;

            int ra = FindRoot(a.WorldIndex);
            int rb = FindRoot(b.WorldIndex);
            if (ra != rb)
                _islandParent[rb] = ra;
        }

        /// <summary>
        /// Groups dynamic bodies connected by contacts or joints. Bodies only fall
        /// asleep or wake up together with the rest of their island.
        /// </summary>
        private void BuildIslands()
        {
            int count = _bodies.Count;
            if (_islandParent.Length < count)
            {
                _islandParent = new int[Math.Max(count, _islandParent.Length * 2)];
                _islandFlag = new byte[_islandParent.Length];
            }

            for (int i = 0; i < count; i++)
                _islandParent[i] = i;

            foreach (var manifold in _manifolds)
            {
                if (manifold.IsActive && !manifold.IsTrigger)
                    Union(manifold.BodyA, manifold.BodyB);
            }

            foreach (var joint in _joints)
            {
                if (joint.IsActive)
                    Union(joint.BodyA, joint.BodyB);
            }
        }

        private bool IsResting(RigidBody2D body)
            => body.Velocity.MagnitudeSquared < SleepVelocityThreshold * SleepVelocityThreshold &&
               Math.Abs(body.AngularVelocity) < SleepAngularVelocityThreshold;

        private static bool IsMovingKinematic(RigidBody2D? body)
            => body != null && body.BodyType == BodyType2D.Kinematic &&
               (body.Velocity.MagnitudeSquared > PhysicsConstants.Epsilon ||
                Math.Abs(body.AngularVelocity) > PhysicsConstants.Epsilon);

        private void FlagIfDrivenByKinematic(RigidBody2D a, RigidBody2D? b, ref bool any)
        {
            if (b == null) return;

            RigidBody2D? driven = null;
            if (IsMovingKinematic(a)) driven = b;
            else if (IsMovingKinematic(b)) driven = a;

            if (driven != null && driven.BodyType == BodyType2D.Dynamic)
            {
                _islandFlag[FindRoot(driven.WorldIndex)] = 1;
                any = true;
            }
        }

        private void WakeDisturbedIslands()
        {
            int count = _bodies.Count;
            Array.Clear(_islandFlag, 0, count);
            bool any = false;

            for (int i = 0; i < count; i++)
            {
                var body = _bodies[i];
                if (body.BodyType != BodyType2D.Dynamic || !body.IsActive || body.IsSleeping)
                    continue;

                if (!IsResting(body))
                {
                    _islandFlag[FindRoot(i)] = 1;
                    any = true;
                }
            }

            foreach (var manifold in _manifolds)
            {
                if (manifold.IsActive && !manifold.IsTrigger)
                    FlagIfDrivenByKinematic(manifold.BodyA, manifold.BodyB, ref any);
            }

            foreach (var joint in _joints)
            {
                if (joint.IsActive)
                    FlagIfDrivenByKinematic(joint.BodyA, joint.BodyB, ref any);
            }

            if (!any) return;

            for (int i = 0; i < count; i++)
            {
                var body = _bodies[i];
                if (body.IsSleeping && body.BodyType == BodyType2D.Dynamic && _islandFlag[FindRoot(i)] == 1)
                    body.WakeUp();
            }
        }

        private void UpdateSleeping(double dt)
        {
            int count = _bodies.Count;

            if (!AllowSleep)
            {
                for (int i = 0; i < count; i++)
                {
                    if (_bodies[i].IsSleeping)
                        _bodies[i].WakeUp();
                }
                return;
            }

            Array.Clear(_islandFlag, 0, count);

            for (int i = 0; i < count; i++)
            {
                var body = _bodies[i];
                if (body.BodyType != BodyType2D.Dynamic || !body.IsActive || body.IsSleeping)
                    continue;

                if (body.CanSleep && IsResting(body))
                    body.SleepTimer += dt;
                else
                    body.SleepTimer = 0;

                if (body.SleepTimer < SleepTimeThreshold)
                    _islandFlag[FindRoot(i)] = 1;
            }

            for (int i = 0; i < count; i++)
            {
                var body = _bodies[i];
                if (body.BodyType == BodyType2D.Dynamic && body.IsActive && !body.IsSleeping &&
                    _islandFlag[FindRoot(i)] == 0)
                {
                    body.Sleep();
                }
            }
        }

        #endregion

        #region Queries

        /// <summary>
        /// Queries all bodies overlapping an AABB (parallel).
        /// </summary>
        public List<RigidBody2D> QueryAABB(AABB2D aabb)
        {
            var results = new ConcurrentBag<RigidBody2D>();

            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (body.IsActive && body.AABB.Overlaps(aabb))
                        results.Add(body);
                });
            }
            finally { _bodiesLock.ExitReadLock(); }

            return results.ToList();
        }

        /// <summary>
        /// Queries all bodies containing a point.
        /// </summary>
        public List<RigidBody2D> QueryPoint(Vector2D point)
        {
            var results = new ConcurrentBag<RigidBody2D>();

            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (!body.IsActive || body.Shape == null) return;
                    if (!body.AABB.Contains(point)) return;

                    if (body.Shape is CircleShape circle)
                    {
                        var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                        if (CollisionDetector2D.PointInCircle(point, center, circle.Radius))
                            results.Add(body);
                    }
                    else if (body.Shape is BoxShape box)
                    {
                        var local = body.WorldToLocal(point);
                        if (Math.Abs(local.X) <= box.HalfWidth && Math.Abs(local.Y) <= box.HalfHeight)
                            results.Add(body);
                    }
                });
            }
            finally { _bodiesLock.ExitReadLock(); }

            return results.ToList();
        }

        /// <summary>
        /// Casts a ray and returns the first hit.
        /// </summary>
        public CollisionDetector2D.RaycastHit2D Raycast(Vector2D origin, Vector2D direction,
            double maxDistance = double.MaxValue, ushort maskBits = 0xFFFF)
        {
            direction = direction.Normalized;
            var closestHit = new CollisionDetector2D.RaycastHit2D { Fraction = double.MaxValue };

            _bodiesLock.EnterReadLock();
            try
            {
                foreach (var body in _bodies)
                {
                    if (!body.IsActive || body.Shape == null) continue;
                    if ((body.CategoryBits & maskBits) == 0) continue;

                    CollisionDetector2D.RaycastHit2D hit;
                    bool didHit = false;

                    if (body.Shape is CircleShape circle)
                    {
                        var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                        didHit = CollisionDetector2D.RaycastCircle(origin, direction, maxDistance,
                            center, circle.Radius, out hit);
                    }
                    else
                    {
                        didHit = CollisionDetector2D.RaycastAABB(origin, direction, maxDistance,
                            body.AABB, out hit);
                    }

                    if (didHit && hit.Fraction < closestHit.Fraction)
                    {
                        closestHit = hit;
                        closestHit.Body = body;
                    }
                }
            }
            finally { _bodiesLock.ExitReadLock(); }

            closestHit.Hit = closestHit.Fraction < double.MaxValue;
            return closestHit;
        }

        /// <summary>
        /// Casts a ray and returns the closest hit.
        /// </summary>
        public RaycastHit2D Raycast(Ray2D ray, ushort maskBits = 0xFFFF)
        {
            return Raycaster2D.Raycast(ray, Bodies, maskBits);
        }

        #endregion

        #region Utility Methods

        /// <summary>
        /// Applies an explosion force at a point (parallel).
        /// </summary>
        public void ApplyExplosionForce(Vector2D center, double radius, double force)
        {
            _bodiesLock.EnterReadLock();
            try
            {
                Parallel.ForEach(_bodies, body =>
                {
                    if (body.BodyType != BodyType2D.Dynamic) return;

                    var delta = body.Position - center;
                    double distance = delta.Magnitude;

                    if (distance > radius || distance < PhysicsConstants.Epsilon) return;

                    double falloff = 1.0 - distance / radius;
                    var direction = delta / distance;
                    body.ApplyImpulse(direction * force * falloff);
                });
            }
            finally { _bodiesLock.ExitReadLock(); }
        }

        /// <summary>
        /// Gets the kinetic energy of all bodies.
        /// </summary>
        public double GetTotalKineticEnergy()
        {
            double total = 0;
            _bodiesLock.EnterReadLock();
            try
            {
                foreach (var body in _bodies)
                {
                    if (body.BodyType == BodyType2D.Dynamic)
                    {
                        total += 0.5 * body.Mass * body.Velocity.MagnitudeSquared;
                        total += 0.5 * body.Inertia * body.AngularVelocity * body.AngularVelocity;
                    }
                }
            }
            finally { _bodiesLock.ExitReadLock(); }
            return total;
        }

        /// <summary>
        /// Resets the simulation.
        /// </summary>
        public void Reset()
        {
            ClearBodies();
            SimulationTime = 0;
        }

        #endregion
    }

    /// <summary>
    /// Continuous collision detection for fast-moving 2D bodies.
    /// </summary>
    public static class ContinuousCollisionDetector2D
    {
        /// <summary>
        /// Detects a collision between two bodies. Kept for API compatibility: tunnelling
        /// is prevented by <see cref="TimeOfImpact"/>, which the world applies while
        /// integrating positions.
        /// </summary>
        public static bool DetectSwept(RigidBody2D bodyA, RigidBody2D bodyB, Manifold2D manifold)
            => CollisionDetector2D.DetectCollision(bodyA, bodyB, manifold);

        /// <summary>
        /// Distance from a point to the surface of a body's shape (0 inside solid shapes).
        /// </summary>
        public static double DistanceToBody(Vector2D point, RigidBody2D body)
        {
            if (body.Shape is CircleShape circle)
            {
                var center = body.Position + Vector2D.Rotate(circle.Offset, body.Rotation);
                return Math.Max(0, (point - center).Magnitude - circle.Radius);
            }

            if (body.Shape is ChainShape chain)
            {
                double best = double.PositiveInfinity;
                for (int i = 0; i < chain.EdgeCount; i++)
                {
                    best = Math.Min(best, DistanceToSegment(point,
                        body.LocalToWorld(chain.Vertices[i]),
                        body.LocalToWorld(chain.Vertices[(i + 1) % chain.VertexCount])));
                }
                return best;
            }

            var vertices = CollisionDetector2D.GetConvexWorldVertices(body);
            if (vertices != null)
            {
                if (vertices.Length >= 3 && CollisionDetector2D.PointInPolygon(point, vertices))
                    return 0;

                double best = double.PositiveInfinity;
                int edges = vertices.Length == 2 ? 1 : vertices.Length;
                for (int i = 0; i < edges; i++)
                    best = Math.Min(best, DistanceToSegment(point, vertices[i], vertices[(i + 1) % vertices.Length]));
                return best;
            }

            var box = body.AABB;
            double dx = Math.Max(0, Math.Max(box.Min.X - point.X, point.X - box.Max.X));
            double dy = Math.Max(0, Math.Max(box.Min.Y - point.Y, point.Y - box.Max.Y));
            return Math.Sqrt(dx * dx + dy * dy);
        }

        private static double DistanceToSegment(Vector2D point, Vector2D a, Vector2D b)
        {
            var ab = b - a;
            double lengthSq = ab.MagnitudeSquared;
            double t = lengthSq > PhysicsConstants.Epsilon ? Math.Clamp(Vector2D.Dot(point - a, ab) / lengthSq, 0, 1) : 0;
            return (point - (a + ab * t)).Magnitude;
        }

        /// <summary>
        /// Conservative advancement of the body's core circle against the other bodies.
        /// Returns the part of <paramref name="dt"/> the body can travel before it would
        /// pass through something; the contact solver handles the impact on the next step.
        /// Bodies that move less than their core radius per step are not swept at all.
        /// </summary>
        /// <param name="body">The moving body (circle or box).</param>
        /// <param name="bodies">All bodies of the world.</param>
        /// <param name="dt">Step duration.</param>
        /// <param name="againstDynamicBodies">Also sweep against other dynamic bodies.</param>
        public static double TimeOfImpact(RigidBody2D body, IReadOnlyList<RigidBody2D> bodies, double dt,
            bool againstDynamicBodies = true)
        {
            double coreRadius;
            var start = body.Position;

            switch (body.Shape)
            {
                case CircleShape circle:
                    coreRadius = circle.Radius;
                    start += Vector2D.Rotate(circle.Offset, body.Rotation);
                    break;
                case BoxShape box:
                    coreRadius = Math.Min(box.HalfWidth, box.HalfHeight);
                    start += Vector2D.Rotate(box.Offset, body.Rotation);
                    break;
                default:
                    return dt;
            }

            var velocity = body.Velocity;
            if (coreRadius <= 0 || velocity.Magnitude * dt <= coreRadius)
                return dt;

            const double tolerance = 1e-4;
            double earliest = dt;

            for (int index = 0; index < bodies.Count; index++)
            {
                var other = bodies[index];
                if (ReferenceEquals(other, body) || !other.IsActive || other.Shape == null || other.IsTrigger)
                    continue;
                if (!againstDynamicBodies && other.BodyType == BodyType2D.Dynamic)
                    continue;
                if (!body.ShouldCollide(other))
                    continue;

                var otherVelocity = other.BodyType == BodyType2D.Static || other.IsSleeping
                    ? Vector2D.Zero
                    : other.Velocity;
                var relativeVelocity = velocity - otherVelocity;
                double relativeSpeed = relativeVelocity.Magnitude;
                if (relativeSpeed < PhysicsConstants.Epsilon)
                    continue;

                // Cheap reject: the swept core must come near the other body's bounds.
                var bounds = other.AABB;
                var end = start + relativeVelocity * dt;
                if (Math.Max(start.X, end.X) + coreRadius < bounds.Min.X || Math.Min(start.X, end.X) - coreRadius > bounds.Max.X ||
                    Math.Max(start.Y, end.Y) + coreRadius < bounds.Min.Y || Math.Min(start.Y, end.Y) - coreRadius > bounds.Max.Y)
                    continue;

                double t = 0;
                for (int iteration = 0; iteration < 32 && t < earliest; iteration++)
                {
                    double distance = DistanceToBody(start + relativeVelocity * t, other) - coreRadius;
                    if (distance <= tolerance)
                    {
                        // Touching at the start of the step is left to the discrete solver.
                        if (iteration > 0)
                        {
                            double overshoot = Math.Min(0.25 * coreRadius, 0.01) / relativeSpeed;
                            earliest = Math.Min(earliest, t + overshoot);
                        }
                        break;
                    }

                    t += distance / relativeSpeed;
                }
            }

            return earliest;
        }
    }
}
