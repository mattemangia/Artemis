using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading;
using Artemis.Collision;
using Artemis.Core;

namespace Artemis.Simulation
{
    /// <summary>
    /// Generic thread-safe object pool for reducing GC pressure in real-time simulations.
    /// </summary>
    /// <typeparam name="T">Type of pooled objects.</typeparam>
    public class ObjectPool<T> where T : class, new()
    {
        private readonly ConcurrentBag<T> _pool;
        private readonly Func<T> _factory;
        private readonly Action<T>? _reset;
        private readonly int _maxSize;
        private int _count;

        /// <summary>
        /// Gets the number of objects currently in the pool.
        /// </summary>
        public int Count => _pool.Count;

        /// <summary>
        /// Gets the total number of objects created by this pool.
        /// </summary>
        public int TotalCreated => _count;

        /// <summary>
        /// Creates a new object pool.
        /// </summary>
        /// <param name="factory">Factory function to create new objects.</param>
        /// <param name="reset">Optional action to reset objects when returned.</param>
        /// <param name="initialSize">Initial pool size.</param>
        /// <param name="maxSize">Maximum pool size.</param>
        public ObjectPool(
            Func<T>? factory = null,
            Action<T>? reset = null,
            int initialSize = 0,
            int maxSize = 10000)
        {
            _pool = new ConcurrentBag<T>();
            _factory = factory ?? (() => new T());
            _reset = reset;
            _maxSize = maxSize;

            // Pre-allocate initial objects
            for (int i = 0; i < initialSize; i++)
            {
                _pool.Add(_factory());
                _count++;
            }
        }

        /// <summary>
        /// Gets an object from the pool or creates a new one.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public T Get()
        {
            if (_pool.TryTake(out T? item))
            {
                return item;
            }

            Interlocked.Increment(ref _count);
            return _factory();
        }

        /// <summary>
        /// Returns an object to the pool.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Return(T item)
        {
            if (_pool.Count < _maxSize)
            {
                _reset?.Invoke(item);
                _pool.Add(item);
            }
        }

        /// <summary>
        /// Clears the pool.
        /// </summary>
        public void Clear()
        {
            while (_pool.TryTake(out _)) { }
        }
    }

    /// <summary>
    /// Specialized pool for lists to avoid allocation.
    /// </summary>
    /// <typeparam name="T">Element type.</typeparam>
    public class ListPool<T>
    {
        private readonly ConcurrentBag<List<T>> _pool;
        private readonly int _initialCapacity;
        private readonly int _maxSize;

        /// <summary>
        /// Creates a new list pool.
        /// </summary>
        public ListPool(int initialCapacity = 16, int maxSize = 1000)
        {
            _pool = new ConcurrentBag<List<T>>();
            _initialCapacity = initialCapacity;
            _maxSize = maxSize;
        }

        /// <summary>
        /// Gets a list from the pool.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public List<T> Get()
        {
            if (_pool.TryTake(out var list))
            {
                return list;
            }
            return new List<T>(_initialCapacity);
        }

        /// <summary>
        /// Returns a list to the pool after clearing it.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void Return(List<T> list)
        {
            if (_pool.Count < _maxSize)
            {
                list.Clear();
                _pool.Add(list);
            }
        }
    }

    /// <summary>
    /// Pool for collision info structures.
    /// </summary>
    public static class CollisionPool
    {
        private static readonly ConcurrentBag<List<CollisionInfo>> _listPool = new();
        private static readonly int MaxPoolSize = 100;

        /// <summary>
        /// Gets a collision list from the pool.
        /// </summary>
        public static List<CollisionInfo> GetList()
        {
            if (_listPool.TryTake(out var list))
            {
                return list;
            }
            return new List<CollisionInfo>(64);
        }

        /// <summary>
        /// Returns a collision list to the pool.
        /// </summary>
        public static void ReturnList(List<CollisionInfo> list)
        {
            if (_listPool.Count < MaxPoolSize)
            {
                list.Clear();
                _listPool.Add(list);
            }
        }
    }

    /// <summary>
    /// Pool for Vector3D arrays (used in contact manifolds, etc.).
    /// </summary>
    public class VectorArrayPool
    {
        private readonly ConcurrentDictionary<int, ConcurrentBag<Vector3D[]>> _pools;
        private readonly int _maxPerSize;

        /// <summary>
        /// Creates a new vector array pool.
        /// </summary>
        public VectorArrayPool(int maxPerSize = 100)
        {
            _pools = new ConcurrentDictionary<int, ConcurrentBag<Vector3D[]>>();
            _maxPerSize = maxPerSize;
        }

        /// <summary>
        /// Gets an array of the specified size.
        /// </summary>
        public Vector3D[] Get(int size)
        {
            var pool = _pools.GetOrAdd(size, _ => new ConcurrentBag<Vector3D[]>());
            if (pool.TryTake(out var array))
            {
                return array;
            }
            return new Vector3D[size];
        }

        /// <summary>
        /// Returns an array to the pool.
        /// </summary>
        public void Return(Vector3D[] array)
        {
            var pool = _pools.GetOrAdd(array.Length, _ => new ConcurrentBag<Vector3D[]>());
            if (pool.Count < _maxPerSize)
            {
                // Clear the array
                Array.Clear(array, 0, array.Length);
                pool.Add(array);
            }
        }
    }

    /// <summary>
    /// Provides pooled temporary allocations for physics calculations.
    /// Uses thread-local storage for lock-free access.
    /// </summary>
    public static class PhysicsAlloc
    {
        [ThreadStatic]
        private static List<Vector3D>? _vectorList;

        [ThreadStatic]
        private static List<double>? _doubleList;

        [ThreadStatic]
        private static List<int>? _intList;

        [ThreadStatic]
        private static Vector3D[]? _vectorBuffer;

        [ThreadStatic]
        private static double[]? _doubleBuffer;

        /// <summary>
        /// Gets a temporary list of vectors (thread-local, cleared on each call).
        /// </summary>
        public static List<Vector3D> TempVectorList
        {
            get
            {
                _vectorList ??= new List<Vector3D>(64);
                _vectorList.Clear();
                return _vectorList;
            }
        }

        /// <summary>
        /// Gets a temporary list of doubles (thread-local, cleared on each call).
        /// </summary>
        public static List<double> TempDoubleList
        {
            get
            {
                _doubleList ??= new List<double>(64);
                _doubleList.Clear();
                return _doubleList;
            }
        }

        /// <summary>
        /// Gets a temporary list of ints (thread-local, cleared on each call).
        /// </summary>
        public static List<int> TempIntList
        {
            get
            {
                _intList ??= new List<int>(64);
                _intList.Clear();
                return _intList;
            }
        }

        /// <summary>
        /// Gets a temporary vector buffer of at least the specified size.
        /// </summary>
        public static Vector3D[] GetVectorBuffer(int minSize)
        {
            if (_vectorBuffer == null || _vectorBuffer.Length < minSize)
            {
                // Round up to power of 2 for efficiency
                int size = 1;
                while (size < minSize) size *= 2;
                _vectorBuffer = new Vector3D[size];
            }
            return _vectorBuffer;
        }

        /// <summary>
        /// Gets a temporary double buffer of at least the specified size.
        /// </summary>
        public static double[] GetDoubleBuffer(int minSize)
        {
            if (_doubleBuffer == null || _doubleBuffer.Length < minSize)
            {
                int size = 1;
                while (size < minSize) size *= 2;
                _doubleBuffer = new double[size];
            }
            return _doubleBuffer;
        }
    }

    /// <summary>
    /// Struct-based temporary scope for pooled objects.
    /// Uses dispose pattern for automatic return.
    /// </summary>
    /// <typeparam name="T">Pooled object type.</typeparam>
    public readonly struct PooledObject<T> : IDisposable where T : class, new()
    {
        private readonly ObjectPool<T> _pool;

        /// <summary>
        /// The pooled object.
        /// </summary>
        public readonly T Value;

        internal PooledObject(ObjectPool<T> pool, T value)
        {
            _pool = pool;
            Value = value;
        }

        /// <summary>
        /// Returns the object to the pool.
        /// </summary>
        public void Dispose()
        {
            _pool.Return(Value);
        }
    }

    /// <summary>
    /// Extensions for object pools.
    /// </summary>
    public static class PoolExtensions
    {
        /// <summary>
        /// Gets a pooled object with automatic disposal.
        /// Usage: using var obj = pool.GetScoped();
        /// </summary>
        public static PooledObject<T> GetScoped<T>(this ObjectPool<T> pool) where T : class, new()
        {
            return new PooledObject<T>(pool, pool.Get());
        }
    }
}
