using System;
using System.Collections.Concurrent;
using System.Linq;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Artemis.Core;
using Artemis.Bodies;
using Artemis.Particles;

namespace Artemis.Simulation
{
    /// <summary>
    /// High-performance precomputation engine for scientific simulations.
    /// Supports multithreading, batch processing, and efficient data storage.
    /// </summary>
    public class PrecomputedSimulation
    {
        #region Fields

        private readonly PhysicsWorld _world;
        private readonly SimulationRecorder _recorder;
        private readonly ConcurrentQueue<Action> _postProcessQueue;
        private volatile bool _isRunning;
        private volatile bool _isPaused;
        private double _progress;

        #endregion

        #region Properties

        /// <summary>
        /// Gets the physics world being simulated.
        /// </summary>
        public PhysicsWorld World => _world;

        /// <summary>
        /// Gets the simulation recorder.
        /// </summary>
        public SimulationRecorder Recorder => _recorder;

        /// <summary>
        /// Gets the simulation progress (0-1).
        /// </summary>
        public double Progress => _progress;

        /// <summary>
        /// Gets whether the simulation is running.
        /// </summary>
        public bool IsRunning => _isRunning;

        /// <summary>
        /// Gets whether the simulation is paused.
        /// </summary>
        public bool IsPaused => _isPaused;

        /// <summary>
        /// Gets or sets the number of substeps per frame.
        /// </summary>
        public int SubSteps { get; set; } = 4;

        /// <summary>
        /// Gets or sets whether to use parallel processing.
        /// </summary>
        public bool UseParallelProcessing { get; set; } = true;

        /// <summary>
        /// Gets or sets the maximum degree of parallelism.
        /// </summary>
        public int MaxParallelism { get; set; } = Environment.ProcessorCount;

        /// <summary>
        /// Event raised on progress update.
        /// </summary>
        public event Action<double>? OnProgress;

        /// <summary>
        /// Event raised when simulation completes.
        /// </summary>
        public event Action<SimulationResult>? OnComplete;

        /// <summary>
        /// Event raised on each frame (for custom processing).
        /// </summary>
        public event Action<int, double>? OnFrame;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a new precomputed simulation.
        /// </summary>
        public PrecomputedSimulation(PhysicsWorld world)
        {
            _world = world;
            _recorder = new SimulationRecorder();
            _postProcessQueue = new ConcurrentQueue<Action>();
            _isRunning = false;
            _isPaused = false;
            _progress = 0;
        }

        /// <summary>
        /// Creates a new precomputed simulation with a fresh world.
        /// </summary>
        public PrecomputedSimulation() : this(new PhysicsWorld())
        {
        }

        #endregion

        #region Precomputation

        /// <summary>
        /// Runs the simulation for a specified duration.
        /// </summary>
        /// <param name="duration">Total simulation time in seconds.</param>
        /// <param name="timeStep">Time step per frame.</param>
        /// <param name="cancellationToken">Cancellation token.</param>
        /// <returns>Simulation result.</returns>
        public async Task<SimulationResult> RunAsync(
            double duration,
            double timeStep = 1.0 / 60.0,
            CancellationToken cancellationToken = default)
        {
            _isRunning = true;
            _isPaused = false;
            _progress = 0;

            var startTime = DateTime.Now;
            int totalFrames = (int)(duration / timeStep);
            int currentFrame = 0;
            double simulatedTime = 0;

            _recorder.Clear();
            _recorder.StartRecording();

            try
            {
                while (simulatedTime < duration && !cancellationToken.IsCancellationRequested)
                {
                    // Handle pause
                    while (_isPaused && !cancellationToken.IsCancellationRequested)
                    {
                        await Task.Delay(100, cancellationToken);
                    }

                    // Simulate frame with substeps
                    double subStepDt = timeStep / SubSteps;
                    for (int sub = 0; sub < SubSteps; sub++)
                    {
                        _world.Update(subStepDt);
                    }

                    // Record frame
                    _recorder.RecordFrame(_world, timeStep);

                    // Fire event
                    OnFrame?.Invoke(currentFrame, simulatedTime);

                    // Process post-processing queue
                    while (_postProcessQueue.TryDequeue(out var action))
                    {
                        action();
                    }

                    simulatedTime += timeStep;
                    currentFrame++;

                    // Update progress
                    _progress = simulatedTime / duration;
                    if (currentFrame % 100 == 0) // Report every 100 frames
                    {
                        OnProgress?.Invoke(_progress);
                    }

                    // Yield occasionally to prevent blocking
                    if (currentFrame % 1000 == 0)
                    {
                        await Task.Yield();
                    }
                }
            }
            finally
            {
                _recorder.StopRecording();
                _isRunning = false;
            }

            var endTime = DateTime.Now;
            var result = new SimulationResult
            {
                Duration = duration,
                TotalFrames = currentFrame,
                TimeStep = timeStep,
                ComputeTime = (endTime - startTime).TotalSeconds,
                Recorder = _recorder,
                EnergyConservationError = _recorder.CalculateEnergyConservationError(),
                AverageFramerate = _recorder.CalculateAverageFramerate(),
                WasCancelled = cancellationToken.IsCancellationRequested
            };

            OnComplete?.Invoke(result);
            return result;
        }

        /// <summary>
        /// Runs the simulation synchronously (blocks).
        /// </summary>
        public SimulationResult Run(double duration, double timeStep = 1.0 / 60.0)
        {
            return RunAsync(duration, timeStep).GetAwaiter().GetResult();
        }

        /// <summary>
        /// Pauses the simulation.
        /// </summary>
        public void Pause()
        {
            _isPaused = true;
        }

        /// <summary>
        /// Resumes the simulation.
        /// </summary>
        public void Resume()
        {
            _isPaused = false;
        }

        /// <summary>
        /// Stops the simulation (cannot be resumed).
        /// </summary>
        public void Stop()
        {
            _isRunning = false;
        }

        #endregion

        #region Batch Processing

        /// <summary>
        /// Runs multiple simulations with parameter variations in parallel.
        /// </summary>
        public async Task<List<SimulationResult>> RunParameterSweepAsync(
            Func<int, PhysicsWorld> worldFactory,
            int variations,
            double duration,
            double timeStep = 1.0 / 60.0,
            CancellationToken cancellationToken = default)
        {
            var results = new ConcurrentBag<SimulationResult>();
            var options = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxParallelism,
                CancellationToken = cancellationToken
            };

            await Task.Run(() =>
            {
                Parallel.For(0, variations, options, i =>
                {
                    var world = worldFactory(i);
                    var sim = new PrecomputedSimulation(world)
                    {
                        SubSteps = SubSteps,
                        UseParallelProcessing = false // Avoid nested parallelism
                    };
                    sim.Recorder.RecordInterval = 10; // Reduce memory for batch

                    var result = sim.Run(duration, timeStep);
                    result.VariationIndex = i;
                    results.Add(result);
                });
            }, cancellationToken);

            return results.OrderBy(r => r.VariationIndex).ToList();
        }

        /// <summary>
        /// Runs a Monte Carlo simulation with random initial conditions.
        /// </summary>
        public async Task<MonteCarloResult> RunMonteCarloAsync(
            Func<Random, PhysicsWorld> worldFactory,
            int samples,
            double duration,
            Func<SimulationResult, double> measurementFunc,
            int? seed = null,
            CancellationToken cancellationToken = default)
        {
            var random = seed.HasValue ? new Random(seed.Value) : new Random();
            var measurements = new ConcurrentBag<double>();
            var options = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxParallelism,
                CancellationToken = cancellationToken
            };

            // Generate seeds for reproducibility
            var seeds = Enumerable.Range(0, samples).Select(_ => random.Next()).ToArray();

            await Task.Run(() =>
            {
                Parallel.For(0, samples, options, i =>
                {
                    var localRandom = new Random(seeds[i]);
                    var world = worldFactory(localRandom);
                    var sim = new PrecomputedSimulation(world)
                    {
                        SubSteps = SubSteps,
                        UseParallelProcessing = false
                    };
                    sim.Recorder.TrackTrajectories = false;
                    sim.Recorder.RecordInterval = 100;

                    var result = sim.Run(duration);
                    double measurement = measurementFunc(result);
                    measurements.Add(measurement);
                });
            }, cancellationToken);

            var values = measurements.ToArray();
            return new MonteCarloResult
            {
                Samples = samples,
                Values = values,
                Mean = values.Average(),
                StandardDeviation = CalculateStandardDeviation(values),
                Min = values.Min(),
                Max = values.Max(),
                Percentile5 = Percentile(values, 5),
                Percentile95 = Percentile(values, 95)
            };
        }

        private static double CalculateStandardDeviation(double[] values)
        {
            if (values.Length < 2) return 0;
            double mean = values.Average();
            double sumSquaredDiff = values.Sum(v => (v - mean) * (v - mean));
            return Math.Sqrt(sumSquaredDiff / (values.Length - 1));
        }

        private static double Percentile(double[] values, double percentile)
        {
            var sorted = values.OrderBy(v => v).ToArray();
            int index = (int)Math.Ceiling(percentile / 100.0 * sorted.Length) - 1;
            return sorted[Math.Max(0, Math.Min(index, sorted.Length - 1))];
        }

        #endregion

        #region Convergence Analysis

        /// <summary>
        /// Runs convergence analysis with different time steps.
        /// </summary>
        public async Task<ConvergenceResult> RunConvergenceAnalysisAsync(
            Func<PhysicsWorld> worldFactory,
            double duration,
            double[] timeSteps,
            Func<SimulationResult, double> errorMetric,
            CancellationToken cancellationToken = default)
        {
            var results = new Dictionary<double, double>();

            foreach (var dt in timeSteps.OrderByDescending(t => t))
            {
                var world = worldFactory();
                var sim = new PrecomputedSimulation(world) { SubSteps = 1 };
                var result = await sim.RunAsync(duration, dt, cancellationToken);
                results[dt] = errorMetric(result);
            }

            // Calculate convergence rate
            var sortedDt = timeSteps.OrderBy(t => t).ToArray();
            double convergenceRate = 0;

            if (sortedDt.Length >= 2)
            {
                double dt1 = sortedDt[0], dt2 = sortedDt[1];
                double e1 = results[dt1], e2 = results[dt2];

                if (e1 > 0 && e2 > 0 && dt1 != dt2)
                {
                    convergenceRate = Math.Log(e2 / e1) / Math.Log(dt2 / dt1);
                }
            }

            return new ConvergenceResult
            {
                TimeSteps = timeSteps,
                Errors = results,
                ConvergenceRate = convergenceRate
            };
        }

        #endregion

        #region Memory-Efficient Processing

        /// <summary>
        /// Runs simulation with streaming output (low memory).
        /// </summary>
        public async Task RunStreamingAsync(
            double duration,
            double timeStep,
            string outputPath,
            CancellationToken cancellationToken = default)
        {
            _isRunning = true;
            int totalFrames = (int)(duration / timeStep);
            double simulatedTime = 0;
            int currentFrame = 0;

            using var stream = new FileStream(outputPath, FileMode.Create, FileAccess.Write, FileShare.None, 65536);
            using var writer = new BinaryWriter(stream);

            // Header
            writer.Write("ARTEMIS_STREAM");
            writer.Write(1); // Version
            writer.Write(totalFrames);
            writer.Write(duration);
            writer.Write(timeStep);
            writer.Write(_world.Bodies.Count);

            // Body IDs
            foreach (var body in _world.Bodies)
            {
                writer.Write(body.Id);
            }

            try
            {
                while (simulatedTime < duration && !cancellationToken.IsCancellationRequested)
                {
                    // Simulate
                    double subStepDt = timeStep / SubSteps;
                    for (int sub = 0; sub < SubSteps; sub++)
                    {
                        _world.Update(subStepDt);
                    }

                    // Write frame directly to disk
                    writer.Write(currentFrame);
                    writer.Write(simulatedTime);

                    foreach (var body in _world.Bodies)
                    {
                        writer.Write(body.Position.X);
                        writer.Write(body.Position.Y);
                        writer.Write(body.Position.Z);
                        writer.Write(body.Velocity.X);
                        writer.Write(body.Velocity.Y);
                        writer.Write(body.Velocity.Z);
                    }

                    simulatedTime += timeStep;
                    currentFrame++;

                    _progress = simulatedTime / duration;
                    if (currentFrame % 1000 == 0)
                    {
                        OnProgress?.Invoke(_progress);
                        await Task.Yield();
                    }
                }
            }
            finally
            {
                _isRunning = false;
            }
        }

        /// <summary>
        /// Reads streaming simulation data frame by frame.
        /// </summary>
        public static IEnumerable<StreamFrame> ReadStreamingData(string filePath)
        {
            using var stream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read, 65536);
            using var reader = new BinaryReader(stream);

            // Header
            string magic = reader.ReadString();
            if (magic != "ARTEMIS_STREAM")
                throw new InvalidDataException("Invalid streaming file");

            int version = reader.ReadInt32();
            int totalFrames = reader.ReadInt32();
            double duration = reader.ReadDouble();
            double timeStep = reader.ReadDouble();
            int bodyCount = reader.ReadInt32();

            var bodyIds = new string[bodyCount];
            for (int i = 0; i < bodyCount; i++)
            {
                bodyIds[i] = reader.ReadString();
            }

            // Frames
            for (int f = 0; f < totalFrames; f++)
            {
                if (stream.Position >= stream.Length)
                    yield break;

                var frame = new StreamFrame
                {
                    FrameIndex = reader.ReadInt32(),
                    Time = reader.ReadDouble(),
                    BodyData = new Dictionary<string, (Vector3D pos, Vector3D vel)>()
                };

                for (int b = 0; b < bodyCount; b++)
                {
                    var pos = new Vector3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
                    var vel = new Vector3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
                    frame.BodyData[bodyIds[b]] = (pos, vel);
                }

                yield return frame;
            }
        }

        #endregion

        #region Queue Post-Processing

        /// <summary>
        /// Queues an action to be executed after the current frame.
        /// </summary>
        public void QueuePostProcess(Action action)
        {
            _postProcessQueue.Enqueue(action);
        }

        #endregion
    }

    #region Result Types

    /// <summary>
    /// Result of a precomputed simulation.
    /// </summary>
    public class SimulationResult
    {
        public double Duration { get; set; }
        public int TotalFrames { get; set; }
        public double TimeStep { get; set; }
        public double ComputeTime { get; set; }
        public SimulationRecorder Recorder { get; set; } = null!;
        public double EnergyConservationError { get; set; }
        public double AverageFramerate { get; set; }
        public bool WasCancelled { get; set; }
        public int VariationIndex { get; set; }

        /// <summary>
        /// Real-time factor (simulated time / compute time).
        /// Values > 1 mean faster than real-time.
        /// </summary>
        public double RealTimeFactor => Duration / ComputeTime;

        /// <summary>
        /// Performance in frames per second.
        /// </summary>
        public double FramesPerSecond => TotalFrames / ComputeTime;
    }

    /// <summary>
    /// Result of Monte Carlo simulation.
    /// </summary>
    public class MonteCarloResult
    {
        public int Samples { get; set; }
        public double[] Values { get; set; } = Array.Empty<double>();
        public double Mean { get; set; }
        public double StandardDeviation { get; set; }
        public double Min { get; set; }
        public double Max { get; set; }
        public double Percentile5 { get; set; }
        public double Percentile95 { get; set; }

        /// <summary>
        /// 95% confidence interval half-width.
        /// </summary>
        public double ConfidenceInterval95 => 1.96 * StandardDeviation / Math.Sqrt(Samples);
    }

    /// <summary>
    /// Result of convergence analysis.
    /// </summary>
    public class ConvergenceResult
    {
        public double[] TimeSteps { get; set; } = Array.Empty<double>();
        public Dictionary<double, double> Errors { get; set; } = new();
        public double ConvergenceRate { get; set; }

        /// <summary>
        /// Gets whether the method shows expected convergence.
        /// </summary>
        public bool IsConverging => ConvergenceRate > 0.5;
    }

    /// <summary>
    /// A frame from streaming data.
    /// </summary>
    public class StreamFrame
    {
        public int FrameIndex { get; set; }
        public double Time { get; set; }
        public Dictionary<string, (Vector3D pos, Vector3D vel)> BodyData { get; set; } = new();
    }

    #endregion
}
