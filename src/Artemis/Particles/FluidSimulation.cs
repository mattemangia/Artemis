using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using System.Numerics;
using Artemis.Bodies;
using Artemis.Core;
using Artemis.Forces;

namespace Artemis.Particles
{
    /// <summary>
    /// High-performance SPH fluid simulation with multi-threading and SIMD.
    /// Simulates realistic fluid behavior with pressure, viscosity, and surface tension.
    /// </summary>
    public class FluidSimulation
    {
        #region Fields

        // Structure of Arrays for SIMD-friendly memory layout
        private float[] _posX;
        private float[] _posY;
        private float[] _posZ;
        private float[] _velX;
        private float[] _velY;
        private float[] _velZ;
        private float[] _forceX;
        private float[] _forceY;
        private float[] _forceZ;
        private float[] _density;
        private float[] _pressure;
        private byte[] _flags; // bit 0 = active

        private int _capacity;
        private int _count;

        // Spatial grid
        private readonly Dictionary<long, List<int>> _grid;
        private readonly float _cellSize;
        private readonly float _invCellSize;

        private readonly List<IPhysicsBody> _colliders;
        private readonly List<IForce> _externalForces;

        // Thread synchronization
        private readonly object _countLock = new object();

        #endregion

        #region Properties

        /// <summary>
        /// Gets all particles as a readonly view.
        /// </summary>
        public IReadOnlyList<FluidParticle> Particles => GetParticleList();

        /// <summary>
        /// Gets the number of particles.
        /// </summary>
        public int ParticleCount => _count;

        /// <summary>
        /// Gets the colliders list.
        /// </summary>
        public IReadOnlyList<IPhysicsBody> Colliders => _colliders;

        /// <summary>
        /// Gets or sets the simulation bounds.
        /// </summary>
        public AABB Bounds { get; set; }

        /// <summary>
        /// Gets or sets the gravity vector.
        /// </summary>
        public Vector3D Gravity { get; set; } = new(0, -9.81, 0);

        /// <summary>
        /// Gets or sets the rest density of the fluid in kg/m³.
        /// </summary>
        public double RestDensity { get; set; } = 1000.0;

        /// <summary>
        /// Gets or sets the gas constant for pressure calculation.
        /// </summary>
        public double GasConstant { get; set; } = 2000.0;

        /// <summary>
        /// Gets or sets the viscosity coefficient.
        /// </summary>
        public double Viscosity { get; set; } = 0.001;

        /// <summary>
        /// Gets or sets the surface tension coefficient.
        /// </summary>
        public double SurfaceTension { get; set; } = 0.0728;

        /// <summary>
        /// Gets or sets the smoothing radius (interaction distance).
        /// </summary>
        public double SmoothingRadius { get; set; } = 0.2;

        /// <summary>
        /// Gets or sets the particle mass.
        /// </summary>
        public double ParticleMass { get; set; } = 0.02;

        /// <summary>
        /// Gets or sets the boundary restitution.
        /// </summary>
        public double BoundaryRestitution { get; set; } = 0.3;

        /// <summary>
        /// Gets or sets the collision restitution with rigid bodies.
        /// </summary>
        public double CollisionRestitution { get; set; } = 0.2;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new fluid simulation.
        /// </summary>
        public FluidSimulation(AABB bounds, double smoothingRadius = 0.2)
        {
            Bounds = bounds;
            SmoothingRadius = smoothingRadius;
            _cellSize = (float)smoothingRadius;
            _invCellSize = 1f / _cellSize;

            _capacity = 1024;
            _count = 0;

            AllocateArrays(_capacity);

            _grid = new Dictionary<long, List<int>>();
            _colliders = new List<IPhysicsBody>();
            _externalForces = new List<IForce>();
        }

        private void AllocateArrays(int capacity)
        {
            _posX = new float[capacity];
            _posY = new float[capacity];
            _posZ = new float[capacity];
            _velX = new float[capacity];
            _velY = new float[capacity];
            _velZ = new float[capacity];
            _forceX = new float[capacity];
            _forceY = new float[capacity];
            _forceZ = new float[capacity];
            _density = new float[capacity];
            _pressure = new float[capacity];
            _flags = new byte[capacity];
        }

        private void EnsureCapacity(int required)
        {
            if (required <= _capacity) return;

            int newCapacity = Math.Max(required, _capacity * 2);

            var newPosX = new float[newCapacity];
            var newPosY = new float[newCapacity];
            var newPosZ = new float[newCapacity];
            var newVelX = new float[newCapacity];
            var newVelY = new float[newCapacity];
            var newVelZ = new float[newCapacity];
            var newForceX = new float[newCapacity];
            var newForceY = new float[newCapacity];
            var newForceZ = new float[newCapacity];
            var newDensity = new float[newCapacity];
            var newPressure = new float[newCapacity];
            var newFlags = new byte[newCapacity];

            Array.Copy(_posX, newPosX, _count);
            Array.Copy(_posY, newPosY, _count);
            Array.Copy(_posZ, newPosZ, _count);
            Array.Copy(_velX, newVelX, _count);
            Array.Copy(_velY, newVelY, _count);
            Array.Copy(_velZ, newVelZ, _count);
            Array.Copy(_forceX, newForceX, _count);
            Array.Copy(_forceY, newForceY, _count);
            Array.Copy(_forceZ, newForceZ, _count);
            Array.Copy(_density, newDensity, _count);
            Array.Copy(_pressure, newPressure, _count);
            Array.Copy(_flags, newFlags, _count);

            _posX = newPosX; _posY = newPosY; _posZ = newPosZ;
            _velX = newVelX; _velY = newVelY; _velZ = newVelZ;
            _forceX = newForceX; _forceY = newForceY; _forceZ = newForceZ;
            _density = newDensity; _pressure = newPressure;
            _flags = newFlags;
            _capacity = newCapacity;
        }

        private List<FluidParticle> GetParticleList()
        {
            var result = new List<FluidParticle>(_count);
            for (int i = 0; i < _count; i++)
            {
                result.Add(new FluidParticle
                {
                    Position = new Vector3D(_posX[i], _posY[i], _posZ[i]),
                    Velocity = new Vector3D(_velX[i], _velY[i], _velZ[i]),
                    Force = new Vector3D(_forceX[i], _forceY[i], _forceZ[i]),
                    Density = _density[i],
                    Pressure = _pressure[i],
                    IsActive = (_flags[i] & 1) != 0
                });
            }
            return result;
        }

        #endregion

        #region Particle Management

        /// <summary>
        /// Adds a fluid particle.
        /// </summary>
        public int AddParticle(Vector3D position, Vector3D velocity = default)
        {
            lock (_countLock)
            {
                EnsureCapacity(_count + 1);

                int index = _count;
                _posX[index] = (float)position.X;
                _posY[index] = (float)position.Y;
                _posZ[index] = (float)position.Z;
                _velX[index] = (float)velocity.X;
                _velY[index] = (float)velocity.Y;
                _velZ[index] = (float)velocity.Z;
                _forceX[index] = 0;
                _forceY[index] = 0;
                _forceZ[index] = 0;
                _density[index] = (float)RestDensity;
                _pressure[index] = 0;
                _flags[index] = 1;

                _count++;
                return index;
            }
        }

        /// <summary>
        /// Adds a block of fluid particles.
        /// </summary>
        public int AddFluidBlock(Vector3D min, Vector3D max, double spacing)
        {
            var positions = new List<(float x, float y, float z)>();
            float halfSpacing = (float)spacing * 0.5f;

            for (float x = (float)min.X + halfSpacing; x < (float)max.X; x += (float)spacing)
            {
                for (float y = (float)min.Y + halfSpacing; y < (float)max.Y; y += (float)spacing)
                {
                    for (float z = (float)min.Z + halfSpacing; z < (float)max.Z; z += (float)spacing)
                    {
                        positions.Add((x, y, z));
                    }
                }
            }

            lock (_countLock)
            {
                EnsureCapacity(_count + positions.Count);

                foreach (var (x, y, z) in positions)
                {
                    int index = _count;
                    _posX[index] = x;
                    _posY[index] = y;
                    _posZ[index] = z;
                    _velX[index] = 0;
                    _velY[index] = 0;
                    _velZ[index] = 0;
                    _forceX[index] = 0;
                    _forceY[index] = 0;
                    _forceZ[index] = 0;
                    _density[index] = (float)RestDensity;
                    _pressure[index] = 0;
                    _flags[index] = 1;
                    _count++;
                }
            }

            return positions.Count;
        }

        /// <summary>
        /// Adds a sphere of fluid particles.
        /// </summary>
        public int AddFluidSphere(Vector3D center, double radius, double spacing)
        {
            var positions = new List<(float x, float y, float z)>();
            float radiusSq = (float)(radius * radius);
            float cx = (float)center.X, cy = (float)center.Y, cz = (float)center.Z;
            float sp = (float)spacing;

            for (float x = cx - (float)radius; x <= cx + (float)radius; x += sp)
            {
                for (float y = cy - (float)radius; y <= cy + (float)radius; y += sp)
                {
                    for (float z = cz - (float)radius; z <= cz + (float)radius; z += sp)
                    {
                        float dx = x - cx, dy = y - cy, dz = z - cz;
                        if (dx * dx + dy * dy + dz * dz <= radiusSq)
                        {
                            positions.Add((x, y, z));
                        }
                    }
                }
            }

            lock (_countLock)
            {
                EnsureCapacity(_count + positions.Count);

                foreach (var (x, y, z) in positions)
                {
                    int index = _count;
                    _posX[index] = x;
                    _posY[index] = y;
                    _posZ[index] = z;
                    _velX[index] = 0;
                    _velY[index] = 0;
                    _velZ[index] = 0;
                    _flags[index] = 1;
                    _density[index] = (float)RestDensity;
                    _count++;
                }
            }

            return positions.Count;
        }

        /// <summary>
        /// Clears all particles.
        /// </summary>
        public void Clear()
        {
            lock (_countLock)
            {
                _count = 0;
                _grid.Clear();
            }
        }

        #endregion

        #region Collider/Force Management

        public void AddCollider(IPhysicsBody body) { if (!_colliders.Contains(body)) _colliders.Add(body); }
        public bool RemoveCollider(IPhysicsBody body) => _colliders.Remove(body);
        public void ClearColliders() => _colliders.Clear();
        public void AddForce(IForce force) { if (!_externalForces.Contains(force)) _externalForces.Add(force); }
        public bool RemoveForce(IForce force) => _externalForces.Remove(force);

        #endregion

        #region Simulation (Multi-threaded)

        /// <summary>
        /// Updates the fluid simulation with multi-threading.
        /// </summary>
        public void Update(double deltaTime, int subSteps = 2)
        {
            float dt = (float)deltaTime / subSteps;

            for (int step = 0; step < subSteps; step++)
            {
                RebuildGrid();
                ComputeDensityPressureParallel();
                ComputeForcesParallel();
                ApplyExternalForcesParallel();
                IntegrateParallel(dt);
                HandleBoundaryCollisionsParallel();
                HandleColliderCollisions();
            }
        }

        private void RebuildGrid()
        {
            _grid.Clear();

            for (int i = 0; i < _count; i++)
            {
                if ((_flags[i] & 1) == 0) continue;

                long cellKey = GetCellKey(_posX[i], _posY[i], _posZ[i]);
                if (!_grid.TryGetValue(cellKey, out var list))
                {
                    list = new List<int>(16);
                    _grid[cellKey] = list;
                }
                list.Add(i);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long GetCellKey(float x, float y, float z)
        {
            int cx = (int)MathF.Floor(x * _invCellSize);
            int cy = (int)MathF.Floor(y * _invCellSize);
            int cz = (int)MathF.Floor(z * _invCellSize);
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (cz & 0x1FFFFF);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private long GetCellKeyFromCoords(int cx, int cy, int cz)
        {
            return ((long)(cx & 0x1FFFFF) << 42) | ((long)(cy & 0x1FFFFF) << 21) | (cz & 0x1FFFFF);
        }

        private void ComputeDensityPressureParallel()
        {
            float h = (float)SmoothingRadius;
            float h2 = h * h;
            float poly6Coeff = (float)(315.0 / (64.0 * Math.PI * Math.Pow(h, 9)));
            float restDensity = (float)RestDensity;
            float gasConstant = (float)GasConstant;
            float particleMass = (float)ParticleMass;

            Parallel.For(0, _count, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                float px = _posX[i], py = _posY[i], pz = _posZ[i];
                float density = 0;

                int cellX = (int)MathF.Floor(px * _invCellSize);
                int cellY = (int)MathF.Floor(py * _invCellSize);
                int cellZ = (int)MathF.Floor(pz * _invCellSize);

                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            long cellKey = GetCellKeyFromCoords(cellX + dx, cellY + dy, cellZ + dz);
                            if (!_grid.TryGetValue(cellKey, out var neighbors)) continue;

                            foreach (int j in neighbors)
                            {
                                if ((_flags[j] & 1) == 0) continue;

                                float ddx = px - _posX[j];
                                float ddy = py - _posY[j];
                                float ddz = pz - _posZ[j];
                                float r2 = ddx * ddx + ddy * ddy + ddz * ddz;

                                if (r2 < h2)
                                {
                                    float diff = h2 - r2;
                                    density += particleMass * poly6Coeff * diff * diff * diff;
                                }
                            }
                        }
                    }
                }

                _density[i] = density;
                _pressure[i] = gasConstant * (density - restDensity);
            });
        }

        private void ComputeForcesParallel()
        {
            float h = (float)SmoothingRadius;
            float h2 = h * h;
            float spikyCoeff = (float)(-45.0 / (Math.PI * Math.Pow(h, 6)));
            float viscosityCoeff = (float)(45.0 / (Math.PI * Math.Pow(h, 6)));
            float particleMass = (float)ParticleMass;
            float viscosity = (float)Viscosity;
            float gravX = (float)Gravity.X;
            float gravY = (float)Gravity.Y;
            float gravZ = (float)Gravity.Z;

            Parallel.For(0, _count, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                float px = _posX[i], py = _posY[i], pz = _posZ[i];
                float vx = _velX[i], vy = _velY[i], vz = _velZ[i];
                float pi_pressure = _pressure[i];
                float pi_density = _density[i];

                float pfX = 0, pfY = 0, pfZ = 0;
                float vfX = 0, vfY = 0, vfZ = 0;

                int cellX = (int)MathF.Floor(px * _invCellSize);
                int cellY = (int)MathF.Floor(py * _invCellSize);
                int cellZ = (int)MathF.Floor(pz * _invCellSize);

                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            long cellKey = GetCellKeyFromCoords(cellX + dx, cellY + dy, cellZ + dz);
                            if (!_grid.TryGetValue(cellKey, out var neighbors)) continue;

                            foreach (int j in neighbors)
                            {
                                if (i == j || (_flags[j] & 1) == 0) continue;

                                float rijX = px - _posX[j];
                                float rijY = py - _posY[j];
                                float rijZ = pz - _posZ[j];
                                float r2 = rijX * rijX + rijY * rijY + rijZ * rijZ;

                                if (r2 < h2 && r2 > 1e-8f)
                                {
                                    float r = MathF.Sqrt(r2);
                                    float invR = 1f / r;
                                    float rijNormX = rijX * invR;
                                    float rijNormY = rijY * invR;
                                    float rijNormZ = rijZ * invR;

                                    float pj_density = _density[j];
                                    float pj_pressure = _pressure[j];

                                    // Pressure force (Spiky gradient)
                                    float diff = h - r;
                                    float pressureGrad = spikyCoeff * diff * diff;
                                    float pressureTerm = particleMass * (pi_pressure + pj_pressure) / (2 * pj_density) * pressureGrad;

                                    pfX -= rijNormX * pressureTerm;
                                    pfY -= rijNormY * pressureTerm;
                                    pfZ -= rijNormZ * pressureTerm;

                                    // Viscosity force (Laplacian)
                                    float viscLaplacian = viscosityCoeff * diff;
                                    float viscTerm = particleMass * viscLaplacian / pj_density * viscosity;

                                    vfX += (_velX[j] - vx) * viscTerm;
                                    vfY += (_velY[j] - vy) * viscTerm;
                                    vfZ += (_velZ[j] - vz) * viscTerm;
                                }
                            }
                        }
                    }
                }

                // Add gravity
                _forceX[i] = pfX + vfX + gravX * pi_density;
                _forceY[i] = pfY + vfY + gravY * pi_density;
                _forceZ[i] = pfZ + vfZ + gravZ * pi_density;
            });
        }

        private void ApplyExternalForcesParallel()
        {
            if (_externalForces.Count == 0) return;

            float particleMass = (float)ParticleMass;

            Parallel.For(0, _count, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                foreach (var force in _externalForces)
                {
                    if (!force.Enabled) continue;

                    var pos = new Vector3D(_posX[i], _posY[i], _posZ[i]);
                    var vel = new Vector3D(_velX[i], _velY[i], _velZ[i]);
                    var f = force.Calculate(pos, vel, particleMass);

                    _forceX[i] += (float)f.X;
                    _forceY[i] += (float)f.Y;
                    _forceZ[i] += (float)f.Z;
                }
            });
        }

        private void IntegrateParallel(float dt)
        {
            int simdWidth = System.Numerics.Vector<float>.Count;
            var dtVec = new System.Numerics.Vector<float>(dt);

            int chunks = (_count + simdWidth - 1) / simdWidth;

            Parallel.For(0, chunks, chunk =>
            {
                int start = chunk * simdWidth;
                int end = Math.Min(start + simdWidth, _count);

                if (end - start == simdWidth && start + simdWidth <= _count)
                {
                    // SIMD path
                    var forceX = new System.Numerics.Vector<float>(_forceX, start);
                    var forceY = new System.Numerics.Vector<float>(_forceY, start);
                    var forceZ = new System.Numerics.Vector<float>(_forceZ, start);
                    var density = new System.Numerics.Vector<float>(_density, start);
                    var velX = new System.Numerics.Vector<float>(_velX, start);
                    var velY = new System.Numerics.Vector<float>(_velY, start);
                    var velZ = new System.Numerics.Vector<float>(_velZ, start);
                    var posX = new System.Numerics.Vector<float>(_posX, start);
                    var posY = new System.Numerics.Vector<float>(_posY, start);
                    var posZ = new System.Numerics.Vector<float>(_posZ, start);

                    // a = F / density
                    var accX = forceX / density;
                    var accY = forceY / density;
                    var accZ = forceZ / density;

                    velX += accX * dtVec;
                    velY += accY * dtVec;
                    velZ += accZ * dtVec;

                    posX += velX * dtVec;
                    posY += velY * dtVec;
                    posZ += velZ * dtVec;

                    velX.CopyTo(_velX, start);
                    velY.CopyTo(_velY, start);
                    velZ.CopyTo(_velZ, start);
                    posX.CopyTo(_posX, start);
                    posY.CopyTo(_posY, start);
                    posZ.CopyTo(_posZ, start);
                }
                else
                {
                    // Scalar fallback
                    for (int i = start; i < end; i++)
                    {
                        if ((_flags[i] & 1) == 0) continue;

                        float invDensity = 1f / _density[i];
                        _velX[i] += _forceX[i] * invDensity * dt;
                        _velY[i] += _forceY[i] * invDensity * dt;
                        _velZ[i] += _forceZ[i] * invDensity * dt;

                        _posX[i] += _velX[i] * dt;
                        _posY[i] += _velY[i] * dt;
                        _posZ[i] += _velZ[i] * dt;
                    }
                }
            });
        }

        private void HandleBoundaryCollisionsParallel()
        {
            float minX = (float)Bounds.Min.X;
            float minY = (float)Bounds.Min.Y;
            float minZ = (float)Bounds.Min.Z;
            float maxX = (float)Bounds.Max.X;
            float maxY = (float)Bounds.Max.Y;
            float maxZ = (float)Bounds.Max.Z;
            float restitution = (float)BoundaryRestitution;
            float radius = (float)SmoothingRadius * 0.5f;

            Parallel.For(0, _count, i =>
            {
                if ((_flags[i] & 1) == 0) return;

                if (_posX[i] - radius < minX)
                {
                    _posX[i] = minX + radius;
                    _velX[i] *= -restitution;
                }
                else if (_posX[i] + radius > maxX)
                {
                    _posX[i] = maxX - radius;
                    _velX[i] *= -restitution;
                }

                if (_posY[i] - radius < minY)
                {
                    _posY[i] = minY + radius;
                    _velY[i] *= -restitution;
                }
                else if (_posY[i] + radius > maxY)
                {
                    _posY[i] = maxY - radius;
                    _velY[i] *= -restitution;
                }

                if (_posZ[i] - radius < minZ)
                {
                    _posZ[i] = minZ + radius;
                    _velZ[i] *= -restitution;
                }
                else if (_posZ[i] + radius > maxZ)
                {
                    _posZ[i] = maxZ - radius;
                    _velZ[i] *= -restitution;
                }
            });
        }

        private void HandleColliderCollisions()
        {
            float radius = (float)SmoothingRadius * 0.5f;
            float restitution = (float)CollisionRestitution;

            foreach (var collider in _colliders)
            {
                if (!collider.IsActive) continue;

                var colliderBounds = collider.BoundingBox.Expanded(radius);

                Parallel.For(0, _count, i =>
                {
                    if ((_flags[i] & 1) == 0) return;

                    var pos = new Vector3D(_posX[i], _posY[i], _posZ[i]);
                    if (!colliderBounds.Contains(pos)) return;

                    if (collider is RigidBody rb)
                    {
                        HandleRigidBodyCollision(i, rb, radius, restitution);
                    }
                });
            }
        }

        private void HandleRigidBodyCollision(int i, RigidBody body, float radius, float restitution)
        {
            var pos = new Vector3D(_posX[i], _posY[i], _posZ[i]);
            var vel = new Vector3D(_velX[i], _velY[i], _velZ[i]);

            switch (body.ShapeType)
            {
                case CollisionShapeType.Sphere:
                    var toParticle = pos - body.Position;
                    double dist = toParticle.Magnitude;
                    if (dist < body.Radius + radius)
                    {
                        var normal = dist > 1e-8 ? toParticle / dist : Vector3D.Up;
                        var closestPoint = body.Position + normal * body.Radius;

                        pos = closestPoint + normal * radius;
                        double velNormal = Vector3D.Dot(vel, normal);
                        if (velNormal < 0)
                        {
                            vel -= normal * velNormal * (1 + restitution);
                        }

                        _posX[i] = (float)pos.X;
                        _posY[i] = (float)pos.Y;
                        _posZ[i] = (float)pos.Z;
                        _velX[i] = (float)vel.X;
                        _velY[i] = (float)vel.Y;
                        _velZ[i] = (float)vel.Z;
                    }
                    break;

                case CollisionShapeType.Box:
                    var localPos = body.Transform.InverseTransformPoint(pos);
                    var halfExtents = body.HalfExtents;

                    var clamped = new Vector3D(
                        Math.Clamp(localPos.X, -halfExtents.X, halfExtents.X),
                        Math.Clamp(localPos.Y, -halfExtents.Y, halfExtents.Y),
                        Math.Clamp(localPos.Z, -halfExtents.Z, halfExtents.Z)
                    );

                    var localDelta = localPos - clamped;
                    double localDist = localDelta.Magnitude;

                    if (localDist < radius || localPos == clamped)
                    {
                        Vector3D normal;
                        if (localDist < 1e-8)
                        {
                            var dists = new double[]
                            {
                                halfExtents.X - Math.Abs(localPos.X),
                                halfExtents.Y - Math.Abs(localPos.Y),
                                halfExtents.Z - Math.Abs(localPos.Z)
                            };

                            int minAxis = 0;
                            for (int a = 1; a < 3; a++)
                                if (dists[a] < dists[minAxis]) minAxis = a;

                            var localNormal = Vector3D.Zero;
                            switch (minAxis)
                            {
                                case 0: localNormal = new Vector3D(Math.Sign(localPos.X), 0, 0); break;
                                case 1: localNormal = new Vector3D(0, Math.Sign(localPos.Y), 0); break;
                                case 2: localNormal = new Vector3D(0, 0, Math.Sign(localPos.Z)); break;
                            }

                            normal = body.Transform.TransformDirection(localNormal);
                            pos = body.Transform.TransformPoint(clamped) + normal * radius;
                        }
                        else
                        {
                            var localNormal = localDelta / localDist;
                            normal = body.Transform.TransformDirection(localNormal);
                            var closestPoint = body.Transform.TransformPoint(clamped);
                            pos = closestPoint + normal * radius;
                        }

                        double velNormal = Vector3D.Dot(vel, normal);
                        if (velNormal < 0)
                        {
                            vel -= normal * velNormal * (1 + restitution);
                        }

                        _posX[i] = (float)pos.X;
                        _posY[i] = (float)pos.Y;
                        _posZ[i] = (float)pos.Z;
                        _velX[i] = (float)vel.X;
                        _velY[i] = (float)vel.Y;
                        _velZ[i] = (float)vel.Z;
                    }
                    break;
            }
        }

        #endregion

        #region Fluent API

        public FluidSimulation WithDensity(double density) { RestDensity = density; return this; }
        public FluidSimulation WithViscosity(double viscosity) { Viscosity = viscosity; return this; }
        public FluidSimulation WithSurfaceTension(double tension) { SurfaceTension = tension; return this; }
        public FluidSimulation WithStiffness(double stiffness) { GasConstant = stiffness; return this; }
        public FluidSimulation WithGravity(Vector3D gravity) { Gravity = gravity; return this; }

        #endregion

        #region Factory Methods

        public static FluidSimulation Water(AABB bounds)
            => new FluidSimulation(bounds).WithDensity(1000).WithViscosity(0.001).WithSurfaceTension(0.0728);

        public static FluidSimulation Oil(AABB bounds)
            => new FluidSimulation(bounds).WithDensity(900).WithViscosity(0.1).WithSurfaceTension(0.03);

        public static FluidSimulation Honey(AABB bounds)
            => new FluidSimulation(bounds).WithDensity(1400).WithViscosity(10.0).WithSurfaceTension(0.05);

        public static FluidSimulation Lava(AABB bounds)
            => new FluidSimulation(bounds).WithDensity(2500).WithViscosity(100.0).WithSurfaceTension(0.4);

        #endregion
    }

    /// <summary>
    /// Represents a fluid particle for SPH simulation.
    /// </summary>
    public struct FluidParticle
    {
        public Vector3D Position;
        public Vector3D Velocity;
        public Vector3D Force;
        public double Density;
        public double Pressure;
        public bool IsActive;
    }
}
