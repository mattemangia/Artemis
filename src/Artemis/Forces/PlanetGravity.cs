using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Represents a gravitational body (planet, moon, asteroid, etc.)
    /// </summary>
    public class GravitationalBody
    {
        public Vector3 Position { get; set; }
        public Vector3 Velocity { get; set; }
        public float Mass { get; set; }
        public float Radius { get; set; }
        public float Density => Mass / (4f / 3f * MathF.PI * Radius * Radius * Radius);
        public bool IsStatic { get; set; }  // If true, doesn't move due to gravity
        public bool HasAtmosphere { get; set; }
        public float AtmosphereHeight { get; set; }
        public float AtmosphereDensity { get; set; }
        public RigidBody? AttachedBody { get; set; }  // Optional rigid body for collisions
        public string Name { get; set; } = "Body";
        public Vector3 Color { get; set; } = Vector3.One;

        /// <summary>
        /// Gravitational parameter (G * M)
        /// </summary>
        public float Mu => PhysicsConstants.G * Mass;

        /// <summary>
        /// Surface gravity
        /// </summary>
        public float SurfaceGravity => Mu / (Radius * Radius);

        /// <summary>
        /// Escape velocity at surface
        /// </summary>
        public float EscapeVelocity => MathF.Sqrt(2 * Mu / Radius);

        /// <summary>
        /// Circular orbital velocity at given altitude
        /// </summary>
        public float OrbitalVelocity(float altitude)
        {
            float r = Radius + altitude;
            return MathF.Sqrt(Mu / r);
        }

        /// <summary>
        /// Orbital period at given altitude (Kepler's third law)
        /// </summary>
        public float OrbitalPeriod(float altitude)
        {
            float r = Radius + altitude;
            return 2 * MathF.PI * MathF.Sqrt(r * r * r / Mu);
        }

        /// <summary>
        /// Sphere of influence radius (Hill sphere approximation)
        /// </summary>
        public float SphereOfInfluence(GravitationalBody parent)
        {
            float distance = (Position - parent.Position).Length();
            return distance * MathF.Pow(Mass / (3 * parent.Mass), 1f / 3f);
        }

        /// <summary>
        /// Calculate atmospheric drag at altitude
        /// </summary>
        public float AtmosphericDrag(float altitude, Vector3 velocity)
        {
            if (!HasAtmosphere || altitude > AtmosphereHeight || altitude < 0)
                return 0;

            float scaleHeight = AtmosphereHeight / 7f;
            float density = AtmosphereDensity * MathF.Exp(-altitude / scaleHeight);
            float speed = velocity.Length();

            // Drag = 0.5 * Cd * rho * v^2 * A (simplified)
            return 0.5f * density * speed * speed;
        }

        // Presets based on real celestial bodies
        public static GravitationalBody Earth() => new()
        {
            Name = "Earth",
            Mass = 5.972e24f,
            Radius = 6.371e6f,
            HasAtmosphere = true,
            AtmosphereHeight = 100000,  // 100 km (Karman line)
            AtmosphereDensity = 1.225f,
            Color = new Vector3(0.2f, 0.4f, 0.8f)
        };

        public static GravitationalBody Moon() => new()
        {
            Name = "Moon",
            Mass = 7.342e22f,
            Radius = 1.737e6f,
            HasAtmosphere = false,
            Color = new Vector3(0.7f, 0.7f, 0.7f)
        };

        public static GravitationalBody Mars() => new()
        {
            Name = "Mars",
            Mass = 6.39e23f,
            Radius = 3.389e6f,
            HasAtmosphere = true,
            AtmosphereHeight = 50000,
            AtmosphereDensity = 0.02f,
            Color = new Vector3(0.8f, 0.3f, 0.2f)
        };

        public static GravitationalBody Sun() => new()
        {
            Name = "Sun",
            Mass = 1.989e30f,
            Radius = 6.96e8f,
            IsStatic = true,
            Color = new Vector3(1f, 0.9f, 0.5f)
        };

        public static GravitationalBody Jupiter() => new()
        {
            Name = "Jupiter",
            Mass = 1.898e27f,
            Radius = 6.991e7f,
            HasAtmosphere = true,
            AtmosphereHeight = 200000,
            AtmosphereDensity = 0.16f,
            Color = new Vector3(0.8f, 0.6f, 0.4f)
        };

        /// <summary>
        /// Create a scaled-down "marble" planet for gameplay
        /// </summary>
        public static GravitationalBody MarblePlanet(float radius, float density = 5500)
        {
            float volume = 4f / 3f * MathF.PI * radius * radius * radius;
            float mass = density * volume;
            return new GravitationalBody
            {
                Name = "Marble",
                Mass = mass,
                Radius = radius
            };
        }

        /// <summary>
        /// Create a game-scale planet (e.g., for games like Mario Galaxy)
        /// </summary>
        public static GravitationalBody GamePlanet(float radius, float surfaceGravity = 9.81f)
        {
            // Calculate mass from desired surface gravity: g = GM/r² => M = gr²/G
            float mass = surfaceGravity * radius * radius / PhysicsConstants.G;
            return new GravitationalBody
            {
                Name = "Planet",
                Mass = mass,
                Radius = radius
            };
        }
    }

    /// <summary>
    /// Manages gravitational interactions between bodies
    /// </summary>
    public class PlanetGravitySystem
    {
        private readonly List<GravitationalBody> _bodies = new();
        private readonly Dictionary<RigidBody, GravitationalBody> _bodyMap = new();

        public float GravitationalConstant { get; set; } = PhysicsConstants.G;
        public bool UseBarnesHut { get; set; } = true;  // Use Barnes-Hut algorithm for large N
        public float BarnesHutTheta { get; set; } = 0.5f;  // Accuracy parameter
        public float SofteningLength { get; set; } = 0.01f;  // Prevents singularity at r=0
        public bool UseParallel { get; set; } = true;
        public int BodyCount => _bodies.Count;

        public IReadOnlyList<GravitationalBody> Bodies => _bodies;

        /// <summary>
        /// Add a gravitational body
        /// </summary>
        public void AddBody(GravitationalBody body)
        {
            _bodies.Add(body);
            if (body.AttachedBody != null)
            {
                _bodyMap[body.AttachedBody] = body;
            }
        }

        /// <summary>
        /// Remove a gravitational body
        /// </summary>
        public void RemoveBody(GravitationalBody body)
        {
            _bodies.Remove(body);
            if (body.AttachedBody != null)
            {
                _bodyMap.Remove(body.AttachedBody);
            }
        }

        /// <summary>
        /// Get gravity at a point from all bodies
        /// </summary>
        public Vector3 GetGravityAt(Vector3 position)
        {
            Vector3 gravity = Vector3.Zero;

            foreach (var body in _bodies)
            {
                Vector3 toBody = body.Position - position;
                float distanceSquared = toBody.LengthSquared() + SofteningLength * SofteningLength;
                float distance = MathF.Sqrt(distanceSquared);

                if (distance > body.Radius * 0.1f)  // Outside the body
                {
                    float gravMag = GravitationalConstant * body.Mass / distanceSquared;
                    gravity += Vector3.Normalize(toBody) * gravMag;
                }
            }

            return gravity;
        }

        /// <summary>
        /// Get the dominant gravitational body at a position (sphere of influence)
        /// </summary>
        public GravitationalBody? GetDominantBody(Vector3 position)
        {
            GravitationalBody? dominant = null;
            float maxInfluence = 0;

            foreach (var body in _bodies)
            {
                float distance = (position - body.Position).Length();
                float influence = body.Mass / (distance * distance);

                if (influence > maxInfluence)
                {
                    maxInfluence = influence;
                    dominant = body;
                }
            }

            return dominant;
        }

        /// <summary>
        /// Update all gravitational bodies (N-body simulation)
        /// </summary>
        public void Update(float deltaTime)
        {
            if (_bodies.Count < 2) return;

            // Calculate accelerations
            var accelerations = new Vector3[_bodies.Count];

            if (UseParallel && _bodies.Count > 100)
            {
                Parallel.For(0, _bodies.Count, i =>
                {
                    accelerations[i] = CalculateAcceleration(_bodies[i]);
                });
            }
            else
            {
                for (int i = 0; i < _bodies.Count; i++)
                {
                    accelerations[i] = CalculateAcceleration(_bodies[i]);
                }
            }

            // Update velocities and positions (leapfrog integration)
            for (int i = 0; i < _bodies.Count; i++)
            {
                if (_bodies[i].IsStatic) continue;

                _bodies[i].Velocity += accelerations[i] * deltaTime;
                _bodies[i].Position += _bodies[i].Velocity * deltaTime;

                // Sync attached rigid body
                if (_bodies[i].AttachedBody != null)
                {
                    _bodies[i].AttachedBody.Position = _bodies[i].Position;
                    _bodies[i].AttachedBody.Velocity = _bodies[i].Velocity;
                }
            }
        }

        private Vector3 CalculateAcceleration(GravitationalBody body)
        {
            Vector3 acceleration = Vector3.Zero;

            foreach (var other in _bodies)
            {
                if (other == body) continue;

                Vector3 toOther = other.Position - body.Position;
                float distanceSquared = toOther.LengthSquared() + SofteningLength * SofteningLength;
                float distance = MathF.Sqrt(distanceSquared);

                float gravMag = GravitationalConstant * other.Mass / distanceSquared;
                acceleration += Vector3.Normalize(toOther) * gravMag;
            }

            return acceleration;
        }

        /// <summary>
        /// Apply gravitational forces to a list of rigid bodies
        /// </summary>
        public void ApplyToRigidBodies(IEnumerable<RigidBody> bodies, float deltaTime)
        {
            foreach (var rb in bodies)
            {
                if (rb.IsStatic) continue;

                Vector3 gravity = GetGravityAt(rb.Position);
                rb.ApplyForce(gravity * rb.Mass);
            }
        }

        /// <summary>
        /// Apply gravitational forces to particles
        /// </summary>
        public void ApplyToParticles(Span<Particles.Particle> particles, float deltaTime)
        {
            foreach (ref var particle in particles)
            {
                if (!particle.IsActive) continue;

                Vector3 gravity = GetGravityAt(particle.Position);
                particle.Velocity += gravity * deltaTime;
            }
        }

        /// <summary>
        /// Check for collisions between gravitational bodies
        /// </summary>
        public List<(GravitationalBody, GravitationalBody)> DetectCollisions()
        {
            var collisions = new List<(GravitationalBody, GravitationalBody)>();

            for (int i = 0; i < _bodies.Count; i++)
            {
                for (int j = i + 1; j < _bodies.Count; j++)
                {
                    float distance = (_bodies[i].Position - _bodies[j].Position).Length();
                    float combinedRadius = _bodies[i].Radius + _bodies[j].Radius;

                    if (distance < combinedRadius)
                    {
                        collisions.Add((_bodies[i], _bodies[j]));
                    }
                }
            }

            return collisions;
        }

        /// <summary>
        /// Merge two colliding bodies (inelastic collision)
        /// </summary>
        public GravitationalBody MergeBodies(GravitationalBody a, GravitationalBody b)
        {
            float totalMass = a.Mass + b.Mass;

            // Conservation of momentum
            Vector3 newVelocity = (a.Mass * a.Velocity + b.Mass * b.Velocity) / totalMass;

            // Center of mass position
            Vector3 newPosition = (a.Mass * a.Position + b.Mass * b.Position) / totalMass;

            // New radius assuming constant density
            float avgDensity = (a.Density + b.Density) / 2;
            float newVolume = totalMass / avgDensity;
            float newRadius = MathF.Pow(3 * newVolume / (4 * MathF.PI), 1f / 3f);

            var merged = new GravitationalBody
            {
                Name = $"{a.Name}+{b.Name}",
                Mass = totalMass,
                Radius = newRadius,
                Position = newPosition,
                Velocity = newVelocity,
                IsStatic = a.IsStatic && b.IsStatic,
                HasAtmosphere = a.HasAtmosphere || b.HasAtmosphere,
                AtmosphereHeight = MathF.Max(a.AtmosphereHeight, b.AtmosphereHeight),
                AtmosphereDensity = (a.AtmosphereDensity + b.AtmosphereDensity) / 2
            };

            RemoveBody(a);
            RemoveBody(b);
            AddBody(merged);

            return merged;
        }

        /// <summary>
        /// Create a stable orbital configuration (e.g., solar system)
        /// </summary>
        public void CreateOrbitalSystem(GravitationalBody central, List<(float distance, float mass, float radius)> orbiters)
        {
            central.IsStatic = true;
            AddBody(central);

            foreach (var (distance, mass, radius) in orbiters)
            {
                var orbiter = new GravitationalBody
                {
                    Mass = mass,
                    Radius = radius,
                    Position = central.Position + new Vector3(distance, 0, 0)
                };

                // Circular orbital velocity
                float orbitalVelocity = MathF.Sqrt(GravitationalConstant * central.Mass / distance);
                orbiter.Velocity = new Vector3(0, 0, orbitalVelocity);

                AddBody(orbiter);
            }
        }

        /// <summary>
        /// Create a planet with satellites
        /// </summary>
        public void CreatePlanetWithMoons(GravitationalBody planet, int moonCount, float minDistance, float maxDistance)
        {
            AddBody(planet);
            var random = new Random();

            for (int i = 0; i < moonCount; i++)
            {
                float distance = minDistance + (float)random.NextDouble() * (maxDistance - minDistance);
                float angle = (float)random.NextDouble() * MathF.PI * 2;
                float inclination = ((float)random.NextDouble() - 0.5f) * 0.2f;  // Slight inclination

                // Moon mass: typically 1/100 to 1/1000 of planet
                float moonMass = planet.Mass * (0.001f + (float)random.NextDouble() * 0.009f);
                float moonRadius = planet.Radius * (0.1f + (float)random.NextDouble() * 0.2f);

                Vector3 position = planet.Position + new Vector3(
                    distance * MathF.Cos(angle) * MathF.Cos(inclination),
                    distance * MathF.Sin(inclination),
                    distance * MathF.Sin(angle) * MathF.Cos(inclination)
                );

                var moon = new GravitationalBody
                {
                    Name = $"Moon {i + 1}",
                    Mass = moonMass,
                    Radius = moonRadius,
                    Position = position
                };

                // Orbital velocity perpendicular to position vector
                float orbitalSpeed = MathF.Sqrt(GravitationalConstant * planet.Mass / distance);
                Vector3 toMoon = moon.Position - planet.Position;
                Vector3 tangent = Vector3.Normalize(Vector3.Cross(toMoon, Vector3.UnitY));
                moon.Velocity = planet.Velocity + tangent * orbitalSpeed;

                AddBody(moon);
            }
        }

        /// <summary>
        /// Calculate Roche limit (distance at which tidal forces exceed self-gravity)
        /// </summary>
        public static float RocheLimit(GravitationalBody primary, GravitationalBody secondary)
        {
            return 2.44f * primary.Radius * MathF.Pow(primary.Density / secondary.Density, 1f / 3f);
        }

        /// <summary>
        /// Calculate tidal force on a body
        /// </summary>
        public Vector3 TidalForce(Vector3 position, float objectRadius, GravitationalBody source)
        {
            Vector3 toSource = source.Position - position;
            float distance = toSource.Length();

            if (distance < source.Radius) return Vector3.Zero;

            // Tidal force gradient
            float tidalMagnitude = 2 * GravitationalConstant * source.Mass * objectRadius /
                                   (distance * distance * distance);

            return Vector3.Normalize(toSource) * tidalMagnitude;
        }
    }

    /// <summary>
    /// Simplified point gravity for game mechanics (e.g., Mario Galaxy style)
    /// </summary>
    public class PointGravity
    {
        public Vector3 Center { get; set; }
        public float Strength { get; set; } = 9.81f;  // Surface gravity
        public float Radius { get; set; } = 1f;       // Planet radius
        public float FalloffRadius { get; set; } = float.MaxValue;  // Gravity cuts off here
        public bool InverseSquare { get; set; } = false;  // If false, constant gravity within range
        public bool PullToSurface { get; set; } = true;   // If true, gravity points to center

        /// <summary>
        /// Get gravity vector at a position
        /// </summary>
        public Vector3 GetGravity(Vector3 position)
        {
            Vector3 toCenter = Center - position;
            float distance = toCenter.Length();

            if (distance > FalloffRadius) return Vector3.Zero;
            if (distance < 0.001f) return Vector3.Zero;

            float gravityMagnitude;
            if (InverseSquare && distance > Radius)
            {
                // Inverse square law outside radius
                float surfaceDistance = Radius;
                gravityMagnitude = Strength * (surfaceDistance * surfaceDistance) / (distance * distance);
            }
            else
            {
                // Constant gravity (game-like)
                gravityMagnitude = Strength;
            }

            return Vector3.Normalize(toCenter) * gravityMagnitude;
        }

        /// <summary>
        /// Apply gravity to a rigid body
        /// </summary>
        public void ApplyTo(RigidBody body)
        {
            Vector3 gravity = GetGravity(body.Position);
            body.ApplyForce(gravity * body.Mass);
        }

        /// <summary>
        /// Check if position is on/near surface
        /// </summary>
        public bool IsOnSurface(Vector3 position, float tolerance = 0.1f)
        {
            float distance = (position - Center).Length();
            return MathF.Abs(distance - Radius) < tolerance;
        }

        /// <summary>
        /// Get "up" direction at a position (opposite of gravity)
        /// </summary>
        public Vector3 GetUpDirection(Vector3 position)
        {
            Vector3 gravity = GetGravity(position);
            return gravity.LengthSquared() > 0 ? -Vector3.Normalize(gravity) : Vector3.UnitY;
        }

        /// <summary>
        /// Project a position onto the surface
        /// </summary>
        public Vector3 ProjectToSurface(Vector3 position)
        {
            Vector3 fromCenter = position - Center;
            if (fromCenter.LengthSquared() < 0.0001f)
            {
                return Center + Vector3.UnitY * Radius;
            }
            return Center + Vector3.Normalize(fromCenter) * Radius;
        }
    }

    /// <summary>
    /// Multi-point gravity for levels with multiple small planets
    /// </summary>
    public class MultiPointGravity
    {
        private readonly List<PointGravity> _gravityPoints = new();

        public IReadOnlyList<PointGravity> GravityPoints => _gravityPoints;

        public void AddPoint(PointGravity point)
        {
            _gravityPoints.Add(point);
        }

        public void RemovePoint(PointGravity point)
        {
            _gravityPoints.Remove(point);
        }

        /// <summary>
        /// Get combined gravity from all points (uses closest planet approach)
        /// </summary>
        public Vector3 GetGravity(Vector3 position, bool useClosest = true)
        {
            if (_gravityPoints.Count == 0) return Vector3.Zero;

            if (useClosest)
            {
                // Find closest gravity point and use only that
                PointGravity? closest = null;
                float minDistance = float.MaxValue;

                foreach (var point in _gravityPoints)
                {
                    float distance = (position - point.Center).Length() - point.Radius;
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        closest = point;
                    }
                }

                return closest?.GetGravity(position) ?? Vector3.Zero;
            }
            else
            {
                // Sum all gravity (more realistic but can be confusing for gameplay)
                Vector3 totalGravity = Vector3.Zero;
                foreach (var point in _gravityPoints)
                {
                    totalGravity += point.GetGravity(position);
                }
                return totalGravity;
            }
        }

        /// <summary>
        /// Get the dominant gravity source
        /// </summary>
        public PointGravity? GetDominantPoint(Vector3 position)
        {
            PointGravity? closest = null;
            float minDistance = float.MaxValue;

            foreach (var point in _gravityPoints)
            {
                float distance = (position - point.Center).Length() - point.Radius;
                if (distance < minDistance)
                {
                    minDistance = distance;
                    closest = point;
                }
            }

            return closest;
        }

        /// <summary>
        /// Apply gravity to rigid bodies
        /// </summary>
        public void ApplyToRigidBodies(IEnumerable<RigidBody> bodies, bool useClosest = true)
        {
            foreach (var body in bodies)
            {
                if (body.IsStatic) continue;
                Vector3 gravity = GetGravity(body.Position, useClosest);
                body.ApplyForce(gravity * body.Mass);
            }
        }
    }
}
