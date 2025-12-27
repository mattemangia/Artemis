using System;
using System.Collections.Generic;
using System.Linq;
using Artemis.Core;

namespace Artemis.Forces
{
    /// <summary>
    /// Force interaction types for combining multiple forces.
    /// </summary>
    public enum ForceBlendMode
    {
        /// <summary>Add all forces together (default physics behavior).</summary>
        Additive,
        /// <summary>Use the strongest force only.</summary>
        Maximum,
        /// <summary>Use the weakest force only.</summary>
        Minimum,
        /// <summary>Average all forces.</summary>
        Average,
        /// <summary>Multiply forces together (for scaling).</summary>
        Multiply,
        /// <summary>Forces cancel each other based on direction.</summary>
        Interference,
        /// <summary>Priority-based: higher priority forces override lower ones.</summary>
        Priority
    }

    /// <summary>
    /// Force interaction rule defining how two forces interact.
    /// </summary>
    public class ForceInteractionRule
    {
        /// <summary>First force type.</summary>
        public Type ForceTypeA { get; set; } = typeof(IForce);

        /// <summary>Second force type.</summary>
        public Type ForceTypeB { get; set; } = typeof(IForce);

        /// <summary>How to blend the forces.</summary>
        public ForceBlendMode BlendMode { get; set; } = ForceBlendMode.Additive;

        /// <summary>Blend factor (0-1) for weighted combinations.</summary>
        public double BlendFactor { get; set; } = 0.5;

        /// <summary>Whether forces amplify each other.</summary>
        public bool Amplify { get; set; }

        /// <summary>Amplification factor when forces align.</summary>
        public double AmplificationFactor { get; set; } = 1.5;

        /// <summary>Whether forces dampen each other.</summary>
        public bool Dampen { get; set; }

        /// <summary>Dampening factor when forces oppose.</summary>
        public double DampeningFactor { get; set; } = 0.5;

        /// <summary>Priority for priority-based blending (higher = more important).</summary>
        public int Priority { get; set; }

        /// <summary>Custom interaction function.</summary>
        public Func<Vector3D, Vector3D, Vector3D>? CustomBlend { get; set; }
    }

    /// <summary>
    /// A composite force that combines multiple forces with interaction rules.
    /// Forces can interact, amplify, dampen, or cancel each other.
    /// </summary>
    public class CompositeForce : IForce
    {
        private readonly List<IForce> _forces = new();
        private readonly List<ForceInteractionRule> _rules = new();
        private readonly Dictionary<(Type, Type), ForceInteractionRule> _ruleCache = new();

        /// <summary>Unique identifier.</summary>
        public string Id { get; set; } = Guid.NewGuid().ToString();

        /// <summary>Whether this force is enabled.</summary>
        public bool Enabled { get; set; } = true;

        /// <summary>Default blend mode when no rule is specified.</summary>
        public ForceBlendMode DefaultBlendMode { get; set; } = ForceBlendMode.Additive;

        /// <summary>Gets the list of forces.</summary>
        public IReadOnlyList<IForce> Forces => _forces;

        /// <summary>Gets the interaction rules.</summary>
        public IReadOnlyList<ForceInteractionRule> Rules => _rules;

        /// <summary>
        /// Adds a force to the composite.
        /// </summary>
        public CompositeForce Add(IForce force)
        {
            _forces.Add(force);
            return this;
        }

        /// <summary>
        /// Adds multiple forces.
        /// </summary>
        public CompositeForce AddRange(IEnumerable<IForce> forces)
        {
            _forces.AddRange(forces);
            return this;
        }

        /// <summary>
        /// Removes a force from the composite.
        /// </summary>
        public bool Remove(IForce force)
        {
            return _forces.Remove(force);
        }

        /// <summary>
        /// Adds an interaction rule.
        /// </summary>
        public CompositeForce AddRule(ForceInteractionRule rule)
        {
            _rules.Add(rule);
            _ruleCache[(rule.ForceTypeA, rule.ForceTypeB)] = rule;
            _ruleCache[(rule.ForceTypeB, rule.ForceTypeA)] = rule;
            return this;
        }

        /// <summary>
        /// Sets a rule for specific force types.
        /// </summary>
        public CompositeForce SetRule<TForceA, TForceB>(ForceBlendMode mode)
            where TForceA : IForce
            where TForceB : IForce
        {
            return AddRule(new ForceInteractionRule
            {
                ForceTypeA = typeof(TForceA),
                ForceTypeB = typeof(TForceB),
                BlendMode = mode
            });
        }

        /// <summary>
        /// Sets forces of the same type to amplify each other.
        /// </summary>
        public CompositeForce SetAmplify<TForce>(double factor = 1.5) where TForce : IForce
        {
            return AddRule(new ForceInteractionRule
            {
                ForceTypeA = typeof(TForce),
                ForceTypeB = typeof(TForce),
                BlendMode = ForceBlendMode.Additive,
                Amplify = true,
                AmplificationFactor = factor
            });
        }

        /// <summary>
        /// Sets opposing forces to dampen each other.
        /// </summary>
        public CompositeForce SetDampen<TForceA, TForceB>(double factor = 0.5)
            where TForceA : IForce
            where TForceB : IForce
        {
            return AddRule(new ForceInteractionRule
            {
                ForceTypeA = typeof(TForceA),
                ForceTypeB = typeof(TForceB),
                Dampen = true,
                DampeningFactor = factor
            });
        }

        /// <summary>
        /// Calculates the combined force with all interactions.
        /// </summary>
        public Vector3D Calculate(Vector3D position, Vector3D velocity, double mass)
        {
            if (!Enabled || _forces.Count == 0)
                return Vector3D.Zero;

            // Calculate individual forces
            var forceVectors = new List<(IForce force, Vector3D vector)>();
            foreach (var force in _forces)
            {
                if (!force.Enabled) continue;
                var f = force.Calculate(position, velocity, mass);
                if (f.LengthSquared > 0)
                {
                    forceVectors.Add((force, f));
                }
            }

            if (forceVectors.Count == 0)
                return Vector3D.Zero;

            if (forceVectors.Count == 1)
                return forceVectors[0].vector;

            // Apply interaction rules
            return CombineWithRules(forceVectors);
        }

        private Vector3D CombineWithRules(List<(IForce force, Vector3D vector)> forces)
        {
            var result = Vector3D.Zero;
            var processed = new HashSet<int>();

            for (int i = 0; i < forces.Count; i++)
            {
                if (processed.Contains(i)) continue;

                var (forceA, vecA) = forces[i];
                var combined = vecA;
                processed.Add(i);

                for (int j = i + 1; j < forces.Count; j++)
                {
                    if (processed.Contains(j)) continue;

                    var (forceB, vecB) = forces[j];
                    var rule = GetRule(forceA.GetType(), forceB.GetType());

                    if (rule != null)
                    {
                        combined = ApplyRule(combined, vecB, rule);
                        processed.Add(j);
                    }
                }

                // Add to result with default blend mode for unmatched forces
                switch (DefaultBlendMode)
                {
                    case ForceBlendMode.Additive:
                        result += combined;
                        break;
                    case ForceBlendMode.Maximum:
                        if (combined.LengthSquared > result.LengthSquared)
                            result = combined;
                        break;
                    case ForceBlendMode.Average:
                        result += combined;
                        break;
                    default:
                        result += combined;
                        break;
                }
            }

            // Handle any remaining unprocessed forces
            for (int i = 0; i < forces.Count; i++)
            {
                if (!processed.Contains(i))
                {
                    result += forces[i].vector;
                }
            }

            if (DefaultBlendMode == ForceBlendMode.Average && processed.Count > 0)
            {
                result /= processed.Count;
            }

            return result;
        }

        private ForceInteractionRule? GetRule(Type typeA, Type typeB)
        {
            if (_ruleCache.TryGetValue((typeA, typeB), out var rule))
                return rule;

            // Check for interface/base class matches
            foreach (var r in _rules)
            {
                if ((r.ForceTypeA.IsAssignableFrom(typeA) && r.ForceTypeB.IsAssignableFrom(typeB)) ||
                    (r.ForceTypeA.IsAssignableFrom(typeB) && r.ForceTypeB.IsAssignableFrom(typeA)))
                {
                    return r;
                }
            }

            return null;
        }

        private Vector3D ApplyRule(Vector3D a, Vector3D b, ForceInteractionRule rule)
        {
            // Custom blend takes priority
            if (rule.CustomBlend != null)
                return rule.CustomBlend(a, b);

            // Check for amplification/dampening based on alignment
            if (rule.Amplify || rule.Dampen)
            {
                double dot = Vector3D.Dot(a.Normalized(), b.Normalized());

                if (rule.Amplify && dot > 0.5) // Aligned forces amplify
                {
                    return (a + b) * rule.AmplificationFactor;
                }
                else if (rule.Dampen && dot < -0.5) // Opposing forces dampen
                {
                    return (a + b) * rule.DampeningFactor;
                }
            }

            // Apply blend mode
            switch (rule.BlendMode)
            {
                case ForceBlendMode.Additive:
                    return a + b;

                case ForceBlendMode.Maximum:
                    return a.LengthSquared > b.LengthSquared ? a : b;

                case ForceBlendMode.Minimum:
                    return a.LengthSquared < b.LengthSquared ? a : b;

                case ForceBlendMode.Average:
                    return (a + b) * 0.5;

                case ForceBlendMode.Multiply:
                    return new Vector3D(a.X * b.X, a.Y * b.Y, a.Z * b.Z);

                case ForceBlendMode.Interference:
                    // Forces that oppose cancel out
                    double alignment = Vector3D.Dot(a.Normalized(), b.Normalized());
                    if (alignment < 0)
                    {
                        double cancellation = Math.Abs(alignment);
                        return (a + b) * (1 - cancellation * 0.8);
                    }
                    return a + b;

                case ForceBlendMode.Priority:
                    // Handled by sorting elsewhere
                    return a + b;

                default:
                    return a + b;
            }
        }
    }

    /// <summary>
    /// System for managing force interactions across an entire simulation.
    /// </summary>
    public class ForceInteractionSystem
    {
        private readonly List<IForce> _globalForces = new();
        private readonly Dictionary<string, CompositeForce> _forceGroups = new();
        private readonly List<ForceInteractionRule> _globalRules = new();
        private readonly Dictionary<string, List<IForce>> _layerForces = new();

        /// <summary>Gets global forces applied to all bodies.</summary>
        public IReadOnlyList<IForce> GlobalForces => _globalForces;

        /// <summary>Gets the force groups.</summary>
        public IReadOnlyDictionary<string, CompositeForce> ForceGroups => _forceGroups;

        /// <summary>
        /// Adds a global force that applies to all bodies.
        /// </summary>
        public ForceInteractionSystem AddGlobalForce(IForce force)
        {
            _globalForces.Add(force);
            return this;
        }

        /// <summary>
        /// Creates a named force group for organizing forces.
        /// </summary>
        public CompositeForce CreateGroup(string name)
        {
            var group = new CompositeForce { Id = name };
            _forceGroups[name] = group;
            return group;
        }

        /// <summary>
        /// Gets or creates a force group.
        /// </summary>
        public CompositeForce GetOrCreateGroup(string name)
        {
            if (!_forceGroups.TryGetValue(name, out var group))
            {
                group = CreateGroup(name);
            }
            return group;
        }

        /// <summary>
        /// Adds a force to a named layer (for layered force systems).
        /// </summary>
        public ForceInteractionSystem AddToLayer(string layer, IForce force)
        {
            if (!_layerForces.TryGetValue(layer, out var forces))
            {
                forces = new List<IForce>();
                _layerForces[layer] = forces;
            }
            forces.Add(force);
            return this;
        }

        /// <summary>
        /// Adds a global interaction rule.
        /// </summary>
        public ForceInteractionSystem AddGlobalRule(ForceInteractionRule rule)
        {
            _globalRules.Add(rule);
            return this;
        }

        /// <summary>
        /// Sets wind and fire forces to interact (wind fans flames).
        /// </summary>
        public ForceInteractionSystem EnableWindFireInteraction(double amplification = 2.0)
        {
            return AddGlobalRule(new ForceInteractionRule
            {
                ForceTypeA = typeof(WindForce),
                ForceTypeB = typeof(BuoyancyForce), // Fire uses buoyancy
                Amplify = true,
                AmplificationFactor = amplification
            });
        }

        /// <summary>
        /// Sets magnetic forces to interfere with each other.
        /// </summary>
        public ForceInteractionSystem EnableMagneticInterference()
        {
            return AddGlobalRule(new ForceInteractionRule
            {
                ForceTypeA = typeof(MagneticForce),
                ForceTypeB = typeof(MagneticForce),
                BlendMode = ForceBlendMode.Interference
            });
        }

        /// <summary>
        /// Sets gravity forces to be additive (multiple gravity sources).
        /// </summary>
        public ForceInteractionSystem EnableMultiGravity()
        {
            return AddGlobalRule(new ForceInteractionRule
            {
                ForceTypeA = typeof(GravityForce),
                ForceTypeB = typeof(PointGravityForce),
                BlendMode = ForceBlendMode.Additive
            });
        }

        /// <summary>
        /// Calculates total force at a position considering all forces and interactions.
        /// </summary>
        public Vector3D CalculateTotalForce(Vector3D position, Vector3D velocity, double mass)
        {
            var allForces = new List<(IForce, Vector3D)>();

            // Collect global forces
            foreach (var force in _globalForces)
            {
                if (force.Enabled)
                {
                    var f = force.Calculate(position, velocity, mass);
                    allForces.Add((force, f));
                }
            }

            // Collect group forces
            foreach (var group in _forceGroups.Values)
            {
                if (group.Enabled)
                {
                    var f = group.Calculate(position, velocity, mass);
                    allForces.Add((group, f));
                }
            }

            // Collect layer forces
            foreach (var layer in _layerForces.Values)
            {
                foreach (var force in layer)
                {
                    if (force.Enabled)
                    {
                        var f = force.Calculate(position, velocity, mass);
                        allForces.Add((force, f));
                    }
                }
            }

            // Apply global rules and sum
            return CombineForces(allForces);
        }

        private Vector3D CombineForces(List<(IForce force, Vector3D vector)> forces)
        {
            if (forces.Count == 0) return Vector3D.Zero;
            if (forces.Count == 1) return forces[0].vector;

            var result = Vector3D.Zero;

            // Simple pairwise interaction
            for (int i = 0; i < forces.Count; i++)
            {
                var (forceA, vecA) = forces[i];
                var contribution = vecA;

                for (int j = i + 1; j < forces.Count; j++)
                {
                    var (forceB, vecB) = forces[j];
                    var rule = FindRule(forceA.GetType(), forceB.GetType());

                    if (rule != null && (rule.Amplify || rule.Dampen))
                    {
                        double dot = Vector3D.Dot(vecA.Normalized(), vecB.Normalized());

                        if (rule.Amplify && dot > 0.5)
                        {
                            contribution *= rule.AmplificationFactor;
                        }
                        else if (rule.Dampen && dot < -0.5)
                        {
                            contribution *= rule.DampeningFactor;
                        }
                    }
                }

                result += contribution;
            }

            return result;
        }

        private ForceInteractionRule? FindRule(Type a, Type b)
        {
            foreach (var rule in _globalRules)
            {
                if ((rule.ForceTypeA.IsAssignableFrom(a) && rule.ForceTypeB.IsAssignableFrom(b)) ||
                    (rule.ForceTypeA.IsAssignableFrom(b) && rule.ForceTypeB.IsAssignableFrom(a)))
                {
                    return rule;
                }
            }
            return null;
        }

        /// <summary>
        /// Creates preset force interaction configurations.
        /// </summary>
        public static class Presets
        {
            /// <summary>
            /// Creates a weather system with interacting wind and forces.
            /// </summary>
            public static ForceInteractionSystem Weather()
            {
                var system = new ForceInteractionSystem();
                system.EnableWindFireInteraction();

                var wind = system.CreateGroup("wind");
                wind.DefaultBlendMode = ForceBlendMode.Additive;

                var thermal = system.CreateGroup("thermal");
                thermal.DefaultBlendMode = ForceBlendMode.Additive;

                return system;
            }

            /// <summary>
            /// Creates a space environment with interacting gravity sources.
            /// </summary>
            public static ForceInteractionSystem Space()
            {
                var system = new ForceInteractionSystem();
                system.EnableMultiGravity();

                var gravity = system.CreateGroup("gravity");
                gravity.DefaultBlendMode = ForceBlendMode.Additive;

                return system;
            }

            /// <summary>
            /// Creates a fluid environment with drag and buoyancy.
            /// </summary>
            public static ForceInteractionSystem Fluid()
            {
                var system = new ForceInteractionSystem();

                system.AddGlobalRule(new ForceInteractionRule
                {
                    ForceTypeA = typeof(DragForce),
                    ForceTypeB = typeof(BuoyancyForce),
                    BlendMode = ForceBlendMode.Additive
                });

                return system;
            }

            /// <summary>
            /// Creates an electromagnetic environment.
            /// </summary>
            public static ForceInteractionSystem Electromagnetic()
            {
                var system = new ForceInteractionSystem();
                system.EnableMagneticInterference();

                var magnetic = system.CreateGroup("magnetic");
                magnetic.SetAmplify<MagneticForce>(1.2);

                return system;
            }
        }
    }
}
