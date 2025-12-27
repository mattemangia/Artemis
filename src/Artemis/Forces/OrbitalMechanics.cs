using System;
using System.Collections.Generic;
using Artemis.Bodies;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Orbital elements for Kepler orbits.
    /// </summary>
    public struct OrbitalElements
    {
        /// <summary>Semi-major axis (average radius).</summary>
        public double SemiMajorAxis;

        /// <summary>Eccentricity (0 = circle, 0-1 = ellipse, 1 = parabola, >1 = hyperbola).</summary>
        public double Eccentricity;

        /// <summary>Inclination in radians.</summary>
        public double Inclination;

        /// <summary>Longitude of ascending node in radians.</summary>
        public double LongitudeAscendingNode;

        /// <summary>Argument of periapsis in radians.</summary>
        public double ArgumentOfPeriapsis;

        /// <summary>Mean anomaly at epoch in radians.</summary>
        public double MeanAnomalyAtEpoch;

        /// <summary>Orbital period in seconds.</summary>
        public double Period;

        /// <summary>Periapsis (closest approach) distance.</summary>
        public double Periapsis => SemiMajorAxis * (1 - Eccentricity);

        /// <summary>Apoapsis (farthest point) distance.</summary>
        public double Apoapsis => SemiMajorAxis * (1 + Eccentricity);

        /// <summary>Whether this is a closed (bound) orbit.</summary>
        public bool IsBound => Eccentricity < 1.0;

        /// <summary>Whether this is a circular orbit.</summary>
        public bool IsCircular => Eccentricity < 0.01;
    }

    /// <summary>
    /// Orbital body representing a celestial object with mass.
    /// </summary>
    public class CelestialBody
    {
        /// <summary>Unique identifier.</summary>
        public string Id { get; set; } = Guid.NewGuid().ToString();

        /// <summary>Display name.</summary>
        public string Name { get; set; } = "Body";

        /// <summary>Position in space.</summary>
        public Vector3D Position { get; set; }

        /// <summary>Velocity.</summary>
        public Vector3D Velocity { get; set; }

        /// <summary>Mass in kg.</summary>
        public double Mass { get; set; }

        /// <summary>Radius in meters.</summary>
        public double Radius { get; set; }

        /// <summary>Atmosphere height in meters (0 if none).</summary>
        public double AtmosphereHeight { get; set; }

        /// <summary>Atmosphere density at the surface.</summary>
        public double AtmosphereDensity { get; set; }

        /// <summary>Visual color (ARGB).</summary>
        public uint Color { get; set; } = 0xFFFFFFFF;

        /// <summary>Parent body (for moons, etc.).</summary>
        public CelestialBody? Parent { get; set; }

        /// <summary>Whether this is a fixed body (like a sun).</summary>
        public bool IsFixed { get; set; }

        /// <summary>Gravitational parameter (G * M).</summary>
        public double Mu => PhysicsConstants.GravitationalConstant * Mass;

        /// <summary>Surface gravity (m/s²).</summary>
        public double SurfaceGravity => Mass > 0 && Radius > 0
            ? PhysicsConstants.GravitationalConstant * Mass / (Radius * Radius)
            : 0;

        /// <summary>Escape velocity at surface.</summary>
        public double EscapeVelocity => Radius > 0
            ? Math.Sqrt(2 * Mu / Radius)
            : 0;
    }

    /// <summary>
    /// Orbital mechanics simulation for space physics.
    /// Supports N-body simulation and Kepler orbit propagation.
    /// </summary>
    public class OrbitalMechanics
    {
        #region Fields

        private readonly List<CelestialBody> _bodies;
        private readonly Dictionary<string, OrbitalElements> _orbits;
        private readonly Dictionary<string, List<Vector3D>> _trajectoryHistory;
        private double _time;

        #endregion

        #region Properties

        /// <summary>Gets all celestial bodies.</summary>
        public IReadOnlyList<CelestialBody> Bodies => _bodies;

        /// <summary>Gets or sets whether to use N-body simulation (true) or Kepler (false).</summary>
        public bool UseNBody { get; set; } = true;

        /// <summary>Gets or sets the softening parameter to prevent singularities.</summary>
        public double Softening { get; set; } = 1e6;

        /// <summary>Gets or sets whether to record trajectory history.</summary>
        public bool RecordTrajectory { get; set; } = false;

        /// <summary>Gets or sets max trajectory history points.</summary>
        public int MaxTrajectoryPoints { get; set; } = 1000;

        /// <summary>Current simulation time.</summary>
        public double Time => _time;

        /// <summary>Time scale multiplier.</summary>
        public double TimeScale { get; set; } = 1.0;

        #endregion

        #region Constructors

        /// <summary>
        /// Creates a new orbital mechanics system.
        /// </summary>
        public OrbitalMechanics()
        {
            _bodies = new List<CelestialBody>();
            _orbits = new Dictionary<string, OrbitalElements>();
            _trajectoryHistory = new Dictionary<string, List<Vector3D>>();
        }

        #endregion

        #region Body Management

        /// <summary>
        /// Adds a celestial body.
        /// </summary>
        public void AddBody(CelestialBody body)
        {
            _bodies.Add(body);
            _trajectoryHistory[body.Id] = new List<Vector3D>();
        }

        /// <summary>
        /// Removes a celestial body.
        /// </summary>
        public bool RemoveBody(CelestialBody body)
        {
            _trajectoryHistory.Remove(body.Id);
            _orbits.Remove(body.Id);
            return _bodies.Remove(body);
        }

        /// <summary>
        /// Gets a body by ID.
        /// </summary>
        public CelestialBody? GetBody(string id)
        {
            return _bodies.Find(b => b.Id == id);
        }

        /// <summary>
        /// Clears all bodies.
        /// </summary>
        public void Clear()
        {
            _bodies.Clear();
            _orbits.Clear();
            _trajectoryHistory.Clear();
            _time = 0;
        }

        #endregion

        #region Update

        /// <summary>
        /// Updates the orbital simulation.
        /// </summary>
        public void Update(double deltaTime)
        {
            deltaTime *= TimeScale;
            if (deltaTime <= 0) return;

            if (UseNBody)
            {
                UpdateNBody(deltaTime);
            }
            else
            {
                UpdateKepler(deltaTime);
            }

            _time += deltaTime;

            // Record trajectories
            if (RecordTrajectory)
            {
                foreach (var body in _bodies)
                {
                    var history = _trajectoryHistory[body.Id];
                    history.Add(body.Position);
                    if (history.Count > MaxTrajectoryPoints)
                    {
                        history.RemoveAt(0);
                    }
                }
            }
        }

        private void UpdateNBody(double deltaTime)
        {
            int n = _bodies.Count;
            var accelerations = new Vector3D[n];

            // Calculate accelerations (O(n²))
            for (int i = 0; i < n; i++)
            {
                if (_bodies[i].IsFixed) continue;

                for (int j = 0; j < n; j++)
                {
                    if (i == j) continue;

                    var r = _bodies[j].Position - _bodies[i].Position;
                    double distSq = r.MagnitudeSquared + Softening * Softening;
                    double dist = Math.Sqrt(distSq);
                    double forceMag = PhysicsConstants.GravitationalConstant * _bodies[j].Mass / distSq;

                    accelerations[i] += r.Normalized * forceMag;
                }
            }

            // Velocity Verlet integration
            for (int i = 0; i < n; i++)
            {
                if (_bodies[i].IsFixed) continue;

                // Update velocity (half step)
                _bodies[i].Velocity += accelerations[i] * (deltaTime * 0.5);

                // Update position
                _bodies[i].Position += _bodies[i].Velocity * deltaTime;
            }

            // Recalculate accelerations at new positions
            var newAccelerations = new Vector3D[n];
            for (int i = 0; i < n; i++)
            {
                if (_bodies[i].IsFixed) continue;

                for (int j = 0; j < n; j++)
                {
                    if (i == j) continue;

                    var r = _bodies[j].Position - _bodies[i].Position;
                    double distSq = r.MagnitudeSquared + Softening * Softening;
                    double forceMag = PhysicsConstants.GravitationalConstant * _bodies[j].Mass / distSq;

                    newAccelerations[i] += r.Normalized * forceMag;
                }
            }

            // Update velocity (second half step)
            for (int i = 0; i < n; i++)
            {
                if (_bodies[i].IsFixed) continue;
                _bodies[i].Velocity += newAccelerations[i] * (deltaTime * 0.5);
            }
        }

        private void UpdateKepler(double deltaTime)
        {
            // For each body with a parent, propagate using Kepler's equations
            foreach (var body in _bodies)
            {
                if (body.IsFixed || body.Parent == null) continue;

                if (!_orbits.TryGetValue(body.Id, out var elements))
                {
                    // Calculate orbital elements from current state
                    elements = CalculateOrbitalElements(body, body.Parent);
                    _orbits[body.Id] = elements;
                }

                // Propagate mean anomaly
                double meanMotion = 2 * Math.PI / elements.Period;
                double meanAnomaly = elements.MeanAnomalyAtEpoch + meanMotion * _time;

                // Solve Kepler's equation for eccentric anomaly
                double eccentricAnomaly = SolveKeplersEquation(meanAnomaly, elements.Eccentricity);

                // Calculate true anomaly
                double trueAnomaly = 2 * Math.Atan2(
                    Math.Sqrt(1 + elements.Eccentricity) * Math.Sin(eccentricAnomaly / 2),
                    Math.Sqrt(1 - elements.Eccentricity) * Math.Cos(eccentricAnomaly / 2)
                );

                // Calculate radius
                double radius = elements.SemiMajorAxis * (1 - elements.Eccentricity * Math.Cos(eccentricAnomaly));

                // Calculate position in orbital plane
                double xOrbit = radius * Math.Cos(trueAnomaly);
                double yOrbit = radius * Math.Sin(trueAnomaly);

                // Rotate to 3D position
                body.Position = body.Parent.Position + RotateOrbitToSpace(
                    xOrbit, yOrbit,
                    elements.LongitudeAscendingNode,
                    elements.Inclination,
                    elements.ArgumentOfPeriapsis
                );

                // Calculate velocity
                double mu = body.Parent.Mu;
                double p = elements.SemiMajorAxis * (1 - elements.Eccentricity * elements.Eccentricity);
                double h = Math.Sqrt(mu * p);

                double vxOrbit = -mu / h * Math.Sin(trueAnomaly);
                double vyOrbit = mu / h * (elements.Eccentricity + Math.Cos(trueAnomaly));

                body.Velocity = body.Parent.Velocity + RotateOrbitToSpace(
                    vxOrbit, vyOrbit,
                    elements.LongitudeAscendingNode,
                    elements.Inclination,
                    elements.ArgumentOfPeriapsis
                );
            }
        }

        private double SolveKeplersEquation(double M, double e, int maxIterations = 30, double tolerance = 1e-10)
        {
            // Newton-Raphson method
            double E = M;
            for (int i = 0; i < maxIterations; i++)
            {
                double dE = (E - e * Math.Sin(E) - M) / (1 - e * Math.Cos(E));
                E -= dE;
                if (Math.Abs(dE) < tolerance) break;
            }
            return E;
        }

        private Vector3D RotateOrbitToSpace(double x, double y, double omega, double i, double w)
        {
            double cosO = Math.Cos(omega);
            double sinO = Math.Sin(omega);
            double cosI = Math.Cos(i);
            double sinI = Math.Sin(i);
            double cosW = Math.Cos(w);
            double sinW = Math.Sin(w);

            double xSpace = (cosO * cosW - sinO * sinW * cosI) * x + (-cosO * sinW - sinO * cosW * cosI) * y;
            double ySpace = (sinO * cosW + cosO * sinW * cosI) * x + (-sinO * sinW + cosO * cosW * cosI) * y;
            double zSpace = (sinW * sinI) * x + (cosW * sinI) * y;

            return new Vector3D(xSpace, ySpace, zSpace);
        }

        #endregion

        #region Orbital Calculations

        /// <summary>
        /// Calculates orbital elements from position and velocity.
        /// </summary>
        public OrbitalElements CalculateOrbitalElements(CelestialBody body, CelestialBody primary)
        {
            var r = body.Position - primary.Position;
            var v = body.Velocity - primary.Velocity;
            double mu = primary.Mu;

            // Specific angular momentum
            var h = Vector3D.Cross(r, v);
            double hMag = h.Magnitude;

            // Eccentricity vector
            var eVec = Vector3D.Cross(v, h) / mu - r.Normalized;
            double e = eVec.Magnitude;

            // Semi-major axis
            double rMag = r.Magnitude;
            double vMag = v.Magnitude;
            double energy = vMag * vMag / 2 - mu / rMag;
            double a = -mu / (2 * energy);

            // Inclination
            double inc = Math.Acos(h.Z / hMag);

            // Node vector
            var n = Vector3D.Cross(Vector3D.Up, h);
            double nMag = n.Magnitude;

            // Longitude of ascending node
            double omega = 0;
            if (nMag > 1e-10)
            {
                omega = Math.Acos(n.X / nMag);
                if (n.Y < 0) omega = 2 * Math.PI - omega;
            }

            // Argument of periapsis
            double w = 0;
            if (nMag > 1e-10 && e > 1e-10)
            {
                w = Math.Acos(Vector3D.Dot(n, eVec) / (nMag * e));
                if (eVec.Z < 0) w = 2 * Math.PI - w;
            }

            // True anomaly
            double trueAnomaly = 0;
            if (e > 1e-10)
            {
                trueAnomaly = Math.Acos(Vector3D.Dot(eVec, r) / (e * rMag));
                if (Vector3D.Dot(r, v) < 0) trueAnomaly = 2 * Math.PI - trueAnomaly;
            }

            // Eccentric anomaly
            double E = 2 * Math.Atan2(
                Math.Sqrt(1 - e) * Math.Sin(trueAnomaly / 2),
                Math.Sqrt(1 + e) * Math.Cos(trueAnomaly / 2)
            );

            // Mean anomaly
            double M = E - e * Math.Sin(E);

            // Period
            double period = 2 * Math.PI * Math.Sqrt(a * a * a / mu);

            return new OrbitalElements
            {
                SemiMajorAxis = a,
                Eccentricity = e,
                Inclination = inc,
                LongitudeAscendingNode = omega,
                ArgumentOfPeriapsis = w,
                MeanAnomalyAtEpoch = M,
                Period = period
            };
        }

        /// <summary>
        /// Creates a circular orbit around a primary body.
        /// </summary>
        public void SetCircularOrbit(CelestialBody body, CelestialBody primary, double radius, double inclination = 0)
        {
            double orbitalSpeed = Math.Sqrt(primary.Mu / radius);

            // Position on orbital plane
            body.Position = primary.Position + new Vector3D(radius, 0, 0);

            // Velocity perpendicular to position, tilted by inclination
            double vx = 0;
            double vy = orbitalSpeed * Math.Cos(inclination);
            double vz = orbitalSpeed * Math.Sin(inclination);

            body.Velocity = primary.Velocity + new Vector3D(vx, vy, vz);
            body.Parent = primary;

            // Update orbital elements
            _orbits[body.Id] = CalculateOrbitalElements(body, primary);
        }

        /// <summary>
        /// Creates an elliptical orbit.
        /// </summary>
        public void SetEllipticalOrbit(
            CelestialBody body,
            CelestialBody primary,
            double periapsis,
            double apoapsis,
            double inclination = 0,
            double argumentOfPeriapsis = 0)
        {
            double a = (periapsis + apoapsis) / 2;
            double e = (apoapsis - periapsis) / (apoapsis + periapsis);

            // Start at periapsis
            double radius = periapsis;
            double orbitalSpeed = Math.Sqrt(primary.Mu * (2 / radius - 1 / a));

            // Position at periapsis
            var posDir = RotateOrbitToSpace(1, 0, 0, inclination, argumentOfPeriapsis);
            body.Position = primary.Position + posDir * radius;

            // Velocity perpendicular to position
            var velDir = RotateOrbitToSpace(0, 1, 0, inclination, argumentOfPeriapsis);
            body.Velocity = primary.Velocity + velDir * orbitalSpeed;

            body.Parent = primary;
            _orbits[body.Id] = CalculateOrbitalElements(body, primary);
        }

        /// <summary>
        /// Calculates the velocity needed for a Hohmann transfer orbit.
        /// </summary>
        public (double deltaV1, double deltaV2) CalculateHohmannTransfer(
            CelestialBody primary,
            double startRadius,
            double endRadius)
        {
            double mu = primary.Mu;

            // Transfer orbit semi-major axis
            double a = (startRadius + endRadius) / 2;

            // Velocities
            double v1 = Math.Sqrt(mu / startRadius);
            double vTransfer1 = Math.Sqrt(mu * (2 / startRadius - 1 / a));
            double deltaV1 = vTransfer1 - v1;

            double v2 = Math.Sqrt(mu / endRadius);
            double vTransfer2 = Math.Sqrt(mu * (2 / endRadius - 1 / a));
            double deltaV2 = v2 - vTransfer2;

            return (deltaV1, deltaV2);
        }

        /// <summary>
        /// Gets the trajectory history for a body.
        /// </summary>
        public IReadOnlyList<Vector3D> GetTrajectory(string bodyId)
        {
            return _trajectoryHistory.TryGetValue(bodyId, out var history)
                ? history
                : Array.Empty<Vector3D>();
        }

        #endregion

        #region Sphere of Influence

        /// <summary>
        /// Calculates the sphere of influence radius.
        /// </summary>
        public double CalculateSphereOfInfluence(CelestialBody body, CelestialBody primary)
        {
            if (body.Parent != primary) return 0;

            var r = body.Position - primary.Position;
            double a = r.Magnitude; // Approximate as current distance

            return a * Math.Pow(body.Mass / primary.Mass, 0.4);
        }

        /// <summary>
        /// Finds the dominant gravitational body at a position.
        /// </summary>
        public CelestialBody? FindDominantBody(Vector3D position)
        {
            CelestialBody? dominant = null;
            double maxAccel = 0;

            foreach (var body in _bodies)
            {
                var r = body.Position - position;
                double distSq = r.MagnitudeSquared;
                if (distSq < 1) continue;

                double accel = body.Mu / distSq;
                if (accel > maxAccel)
                {
                    maxAccel = accel;
                    dominant = body;
                }
            }

            return dominant;
        }

        #endregion

        #region Presets

        /// <summary>Creates the Sun.</summary>
        public static CelestialBody CreateSun() => new()
        {
            Name = "Sun",
            Mass = 1.989e30,
            Radius = 6.96e8,
            Color = 0xFFFFFF00,
            IsFixed = true
        };

        /// <summary>Creates Earth.</summary>
        public static CelestialBody CreateEarth() => new()
        {
            Name = "Earth",
            Mass = 5.972e24,
            Radius = 6.371e6,
            Color = 0xFF4488FF
        };

        /// <summary>Creates the Moon.</summary>
        public static CelestialBody CreateMoon() => new()
        {
            Name = "Moon",
            Mass = 7.342e22,
            Radius = 1.737e6,
            Color = 0xFFCCCCCC
        };

        /// <summary>Creates Mars.</summary>
        public static CelestialBody CreateMars() => new()
        {
            Name = "Mars",
            Mass = 6.39e23,
            Radius = 3.39e6,
            Color = 0xFFFF6644
        };

        /// <summary>Creates Jupiter.</summary>
        public static CelestialBody CreateJupiter() => new()
        {
            Name = "Jupiter",
            Mass = 1.898e27,
            Radius = 6.991e7,
            Color = 0xFFDDCC88
        };

        /// <summary>Creates a small asteroid.</summary>
        public static CelestialBody CreateAsteroid(double mass = 1e10, double radius = 100) => new()
        {
            Name = "Asteroid",
            Mass = mass,
            Radius = radius,
            Color = 0xFF888888
        };

        #endregion
    }
}
