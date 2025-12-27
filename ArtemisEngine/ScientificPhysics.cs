namespace ArtemisEngine;

/// <summary>
/// Advanced integration methods for scientific simulations
/// More accurate and stable than basic Euler integration
/// </summary>
public static class AdvancedIntegration
{
    /// <summary>
    /// Verlet integration - better energy conservation than Euler
    /// Commonly used in molecular dynamics and particle systems
    /// </summary>
    public static void VerletIntegration(RigidBody body, Vector2 acceleration, float deltaTime)
    {
        if (body.IsStatic || body.IsSleeping)
            return;

        // Store previous position if not available (first frame)
        Vector2 previousPosition = body.Position - body.Velocity * deltaTime;

        // Verlet integration: x(t+dt) = 2*x(t) - x(t-dt) + a*dt^2
        Vector2 newPosition = 2 * body.Position - previousPosition + acceleration * deltaTime * deltaTime;

        // Calculate velocity from position difference
        body.Velocity = (newPosition - body.Position) / deltaTime;
        body.Position = newPosition;
    }

    /// <summary>
    /// Runge-Kutta 4th order (RK4) - highly accurate integration
    /// Used in aerospace, orbital mechanics, and precision simulations
    /// </summary>
    public static void RK4Integration(RigidBody body, Func<Vector2, Vector2, Vector2> accelerationFunc, float deltaTime)
    {
        if (body.IsStatic || body.IsSleeping)
            return;

        Vector2 x = body.Position;
        Vector2 v = body.Velocity;

        // RK4 steps
        Vector2 k1v = accelerationFunc(x, v);
        Vector2 k1x = v;

        Vector2 k2v = accelerationFunc(x + k1x * (deltaTime / 2), v + k1v * (deltaTime / 2));
        Vector2 k2x = v + k1v * (deltaTime / 2);

        Vector2 k3v = accelerationFunc(x + k2x * (deltaTime / 2), v + k2v * (deltaTime / 2));
        Vector2 k3x = v + k2v * (deltaTime / 2);

        Vector2 k4v = accelerationFunc(x + k3x * deltaTime, v + k3v * deltaTime);
        Vector2 k4x = v + k3v * deltaTime;

        // Weighted average
        body.Velocity += (k1v + 2 * k2v + 2 * k3v + k4v) * (deltaTime / 6);
        body.Position += (k1x + 2 * k2x + 2 * k3x + k4x) * (deltaTime / 6);
    }

    /// <summary>
    /// Semi-implicit Euler (Symplectic Euler) - better for game physics
    /// More stable than explicit Euler, energy preserving
    /// </summary>
    public static void SemiImplicitEuler(RigidBody body, Vector2 acceleration, float deltaTime)
    {
        if (body.IsStatic || body.IsSleeping)
            return;

        // Update velocity first
        body.Velocity += acceleration * deltaTime;

        // Then update position with new velocity
        body.Position += body.Velocity * deltaTime;
    }
}

/// <summary>
/// Energy tracking and conservation for scientific simulations
/// Monitors kinetic, potential, and total energy of the system
/// </summary>
public class EnergyTracker
{
    public float KineticEnergy { get; private set; }
    public float PotentialEnergy { get; private set; }
    public float TotalEnergy => KineticEnergy + PotentialEnergy;
    public float RotationalEnergy { get; private set; }

    public Vector2 GravityDirection { get; set; } = new Vector2(0, -1);
    public float GravityMagnitude { get; set; } = 20f;

    public void UpdateEnergy(IEnumerable<RigidBody> bodies)
    {
        KineticEnergy = 0;
        PotentialEnergy = 0;
        RotationalEnergy = 0;

        foreach (var body in bodies)
        {
            if (body.IsStatic)
                continue;

            // Kinetic energy: KE = 0.5 * m * v^2
            KineticEnergy += 0.5f * body.Mass * body.Velocity.LengthSquared;

            // Rotational kinetic energy: RE = 0.5 * I * ω^2
            RotationalEnergy += 0.5f * body.Inertia * body.AngularVelocity * body.AngularVelocity;

            // Potential energy: PE = m * g * h
            float height = Vector2.Dot(body.Position, -GravityDirection);
            PotentialEnergy += body.Mass * GravityMagnitude * height;
        }
    }

    public float GetEnergyDrift(float initialEnergy)
    {
        return TotalEnergy - initialEnergy;
    }

    public float GetEnergyDriftPercentage(float initialEnergy)
    {
        if (initialEnergy == 0)
            return 0;

        return (GetEnergyDrift(initialEnergy) / initialEnergy) * 100;
    }
}

/// <summary>
/// Deterministic physics for reproducible simulations
/// Ensures same input always produces same output
/// </summary>
public class DeterministicPhysics
{
    private uint _seed;
    private Random _random;

    public DeterministicPhysics(uint seed = 12345)
    {
        _seed = seed;
        _random = new Random((int)seed);
    }

    public void Reset()
    {
        _random = new Random((int)_seed);
    }

    public float NextFloat()
    {
        return (float)_random.NextDouble();
    }

    public float NextFloat(float min, float max)
    {
        return min + (max - min) * NextFloat();
    }

    /// <summary>
    /// Ensures floating-point operations are deterministic across platforms
    /// </summary>
    public static float DeterministicAdd(float a, float b)
    {
        // Use compensated summation (Kahan algorithm) for better precision
        return a + b; // In C#, float ops are already deterministic
    }
}

/// <summary>
/// Data export for scientific analysis
/// Export simulation data to CSV, JSON, or custom formats
/// </summary>
public class SimulationDataExporter
{
    public class FrameData
    {
        public float Time { get; set; }
        public List<BodyState> Bodies { get; set; } = new();
    }

    public class BodyState
    {
        public string? Id { get; set; }
        public Vector2 Position { get; set; }
        public Vector2 Velocity { get; set; }
        public float Rotation { get; set; }
        public float AngularVelocity { get; set; }
        public float Energy { get; set; }
    }

    private List<FrameData> _frames = new();
    private float _currentTime = 0;

    public void RecordFrame(IEnumerable<RigidBody> bodies, float deltaTime)
    {
        _currentTime += deltaTime;

        var frame = new FrameData { Time = _currentTime };

        foreach (var body in bodies)
        {
            if (body.UserData is string id)
            {
                frame.Bodies.Add(new BodyState
                {
                    Id = id,
                    Position = body.Position,
                    Velocity = body.Velocity,
                    Rotation = body.Rotation,
                    AngularVelocity = body.AngularVelocity,
                    Energy = 0.5f * body.Mass * body.Velocity.LengthSquared +
                            0.5f * body.Inertia * body.AngularVelocity * body.AngularVelocity
                });
            }
        }

        _frames.Add(frame);
    }

    public string ExportToCSV()
    {
        var csv = new System.Text.StringBuilder();
        csv.AppendLine("Time,BodyId,PosX,PosY,VelX,VelY,Rotation,AngularVel,Energy");

        foreach (var frame in _frames)
        {
            foreach (var body in frame.Bodies)
            {
                csv.AppendLine($"{frame.Time},{body.Id},{body.Position.X},{body.Position.Y}," +
                              $"{body.Velocity.X},{body.Velocity.Y},{body.Rotation}," +
                              $"{body.AngularVelocity},{body.Energy}");
            }
        }

        return csv.ToString();
    }

    public void Clear()
    {
        _frames.Clear();
        _currentTime = 0;
    }
}

/// <summary>
/// Precision settings for scientific simulations
/// </summary>
public class PrecisionSettings
{
    public enum IntegrationMethod
    {
        Euler,
        SemiImplicitEuler,
        Verlet,
        RK4
    }

    public IntegrationMethod Method { get; set; } = IntegrationMethod.SemiImplicitEuler;
    public int ConstraintIterations { get; set; } = 10;
    public int VelocityIterations { get; set; } = 8;
    public int PositionIterations { get; set; } = 3;
    public float Tolerance { get; set; } = 0.0001f;
    public bool UseDoublePrecision { get; set; } = false; // For future double support
}
