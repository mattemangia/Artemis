using System;
using System.Runtime.CompilerServices;
using Artemis.Core;

namespace Artemis.Physics2D
{
    /// <summary>
    /// Base class for 2D physics joints/constraints.
    /// </summary>
    public abstract class Joint2D
    {
        /// <summary>Unique identifier for this joint.</summary>
        public string Id { get; }

        /// <summary>First connected body.</summary>
        public RigidBody2D BodyA { get; protected set; } = null!;

        /// <summary>Second connected body (null for world-anchored joints).</summary>
        public RigidBody2D? BodyB { get; protected set; }

        /// <summary>Whether the joint is active.</summary>
        public bool IsActive { get; set; } = true;

        /// <summary>Whether connected bodies can collide with each other.</summary>
        public bool CollideConnected { get; set; }

        /// <summary>User-defined data.</summary>
        public object? UserData { get; set; }

        protected Joint2D(string? id = null)
        {
            Id = id ?? $"Joint2D_{Guid.NewGuid():N}";
        }

        /// <summary>
        /// Prepares the joint for solving.
        /// </summary>
        public abstract void PreSolve(double dt);

        /// <summary>
        /// Solves the velocity constraint.
        /// </summary>
        public abstract void SolveVelocity();

        /// <summary>
        /// Solves the position constraint.
        /// </summary>
        public abstract bool SolvePosition();

        /// <summary>
        /// Gets the current reaction force.
        /// </summary>
        public abstract Vector2D GetReactionForce(double invDt);

        /// <summary>
        /// Gets the current reaction torque.
        /// </summary>
        public abstract double GetReactionTorque(double invDt);
    }

    /// <summary>
    /// Distance joint - maintains a fixed distance between two anchor points.
    /// </summary>
    public class DistanceJoint2D : Joint2D
    {
        /// <summary>Anchor point on body A in local coordinates.</summary>
        public Vector2D LocalAnchorA { get; set; }

        /// <summary>Anchor point on body B in local coordinates.</summary>
        public Vector2D LocalAnchorB { get; set; }

        /// <summary>Target distance between anchors.</summary>
        public double Length { get; set; }

        /// <summary>Frequency for soft constraint (0 = rigid).</summary>
        public double Frequency { get; set; }

        /// <summary>Damping ratio for soft constraint (0-1).</summary>
        public double DampingRatio { get; set; } = 0.7;

        /// <summary>Minimum distance (0 = no limit).</summary>
        public double MinLength { get; set; }

        /// <summary>Maximum distance (0 = no limit).</summary>
        public double MaxLength { get; set; }

        // Solver temp
        private Vector2D _u;
        private Vector2D _rA;
        private Vector2D _rB;
        private double _mass;
        private double _impulse;
        private double _gamma;
        private double _bias;

        public DistanceJoint2D(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchorA, Vector2D anchorB, string? id = null) : base(id)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = bodyA.WorldToLocal(anchorA);
            LocalAnchorB = bodyB.WorldToLocal(anchorB);
            Length = Vector2D.Distance(anchorA, anchorB);
        }

        public override void PreSolve(double dt)
        {
            if (BodyB == null) return;

            _rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            _rB = Vector2D.Rotate(LocalAnchorB, BodyB.Rotation);

            _u = BodyB.Position + _rB - BodyA.Position - _rA;

            double length = _u.Magnitude;
            if (length > PhysicsConstants.Epsilon)
                _u /= length;
            else
                _u = Vector2D.Zero;

            double crA = Vector2D.Cross(_rA, _u);
            double crB = Vector2D.Cross(_rB, _u);

            double invMass = BodyA.InverseMass + BodyB.InverseMass +
                             BodyA.InverseInertia * crA * crA +
                             BodyB.InverseInertia * crB * crB;

            _mass = invMass > 0 ? 1.0 / invMass : 0;

            if (Frequency > 0)
            {
                double omega = 2 * Math.PI * Frequency;
                double d = 2 * _mass * DampingRatio * omega;
                double k = _mass * omega * omega;

                _gamma = dt * (d + dt * k);
                _gamma = _gamma > 0 ? 1.0 / _gamma : 0;
                _bias = (length - Length) * dt * k * _gamma;
                _mass = invMass + _gamma;
                _mass = _mass > 0 ? 1.0 / _mass : 0;
            }
            else
            {
                _gamma = 0;
                _bias = 0;
            }

            // Warm starting
            var P = _u * _impulse;
            BodyA.Velocity -= P * BodyA.InverseMass;
            BodyA.AngularVelocity -= BodyA.InverseInertia * Vector2D.Cross(_rA, P);
            BodyB.Velocity += P * BodyB.InverseMass;
            BodyB.AngularVelocity += BodyB.InverseInertia * Vector2D.Cross(_rB, P);
        }

        public override void SolveVelocity()
        {
            if (BodyB == null) return;

            var vA = BodyA.Velocity + Vector2D.Cross(BodyA.AngularVelocity, _rA);
            var vB = BodyB.Velocity + Vector2D.Cross(BodyB.AngularVelocity, _rB);

            double Cdot = Vector2D.Dot(_u, vB - vA);
            double impulse = -_mass * (Cdot + _bias + _gamma * _impulse);
            _impulse += impulse;

            var P = _u * impulse;
            BodyA.Velocity -= P * BodyA.InverseMass;
            BodyA.AngularVelocity -= BodyA.InverseInertia * Vector2D.Cross(_rA, P);
            BodyB.Velocity += P * BodyB.InverseMass;
            BodyB.AngularVelocity += BodyB.InverseInertia * Vector2D.Cross(_rB, P);
        }

        public override bool SolvePosition()
        {
            if (Frequency > 0 || BodyB == null)
                return true;

            var rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            var rB = Vector2D.Rotate(LocalAnchorB, BodyB.Rotation);

            var u = BodyB.Position + rB - BodyA.Position - rA;
            double length = u.Magnitude;
            u = length > PhysicsConstants.Epsilon ? u / length : Vector2D.Zero;

            double C = length - Length;
            C = Math.Clamp(C, -0.2, 0.2); // Max correction

            double crA = Vector2D.Cross(rA, u);
            double crB = Vector2D.Cross(rB, u);
            double invMass = BodyA.InverseMass + BodyB.InverseMass +
                             BodyA.InverseInertia * crA * crA +
                             BodyB.InverseInertia * crB * crB;

            double impulse = invMass > 0 ? -C / invMass : 0;
            var P = u * impulse;

            BodyA.Position -= P * BodyA.InverseMass;
            BodyA.Rotation -= BodyA.InverseInertia * Vector2D.Cross(rA, P);
            BodyB.Position += P * BodyB.InverseMass;
            BodyB.Rotation += BodyB.InverseInertia * Vector2D.Cross(rB, P);

            return Math.Abs(C) < 0.005;
        }

        public override Vector2D GetReactionForce(double invDt)
            => _u * (_impulse * invDt);

        public override double GetReactionTorque(double invDt) => 0;
    }

    /// <summary>
    /// Revolute joint - allows rotation around a single point.
    /// </summary>
    public class RevoluteJoint2D : Joint2D
    {
        /// <summary>Anchor point on body A in local coordinates.</summary>
        public Vector2D LocalAnchorA { get; set; }

        /// <summary>Anchor point on body B in local coordinates.</summary>
        public Vector2D LocalAnchorB { get; set; }

        /// <summary>Reference angle between bodies.</summary>
        public double ReferenceAngle { get; set; }

        /// <summary>Whether angle limits are enabled.</summary>
        public bool EnableLimit { get; set; }

        /// <summary>Lower angle limit in radians.</summary>
        public double LowerAngle { get; set; }

        /// <summary>Upper angle limit in radians.</summary>
        public double UpperAngle { get; set; }

        /// <summary>Whether the motor is enabled.</summary>
        public bool EnableMotor { get; set; }

        /// <summary>Motor target speed in radians/second.</summary>
        public double MotorSpeed { get; set; }

        /// <summary>Maximum motor torque.</summary>
        public double MaxMotorTorque { get; set; }

        // Solver temp
        private Vector2D _rA;
        private Vector2D _rB;
        private double _motorMass;
        private double _motorImpulse;
        private Vector2D _impulse;
        private double _limitImpulse;
        private double _limitState;
        private double[,] _mass = new double[2, 2];

        public RevoluteJoint2D(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchor, string? id = null) : base(id)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = bodyA.WorldToLocal(anchor);
            LocalAnchorB = bodyB?.WorldToLocal(anchor) ?? Vector2D.Zero;
            ReferenceAngle = (BodyB?.Rotation ?? 0) - bodyA.Rotation;
        }

        public override void PreSolve(double dt)
        {
            _rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            _rB = BodyB != null ? Vector2D.Rotate(LocalAnchorB, BodyB.Rotation) : Vector2D.Zero;

            double mA = BodyA.InverseMass;
            double mB = BodyB?.InverseMass ?? 0;
            double iA = BodyA.InverseInertia;
            double iB = BodyB?.InverseInertia ?? 0;

            // Mass matrix
            _mass[0, 0] = mA + mB + _rA.Y * _rA.Y * iA + _rB.Y * _rB.Y * iB;
            _mass[0, 1] = -_rA.Y * _rA.X * iA - _rB.Y * _rB.X * iB;
            _mass[1, 0] = _mass[0, 1];
            _mass[1, 1] = mA + mB + _rA.X * _rA.X * iA + _rB.X * _rB.X * iB;

            _motorMass = iA + iB;
            if (_motorMass > 0)
                _motorMass = 1.0 / _motorMass;

            if (!EnableMotor)
                _motorImpulse = 0;

            if (EnableLimit)
            {
                double angle = (BodyB?.Rotation ?? 0) - BodyA.Rotation - ReferenceAngle;
                if (Math.Abs(UpperAngle - LowerAngle) < 2 * PhysicsConstants.Epsilon)
                    _limitState = 0; // Equal limits
                else if (angle <= LowerAngle)
                    _limitState = -1; // At lower limit
                else if (angle >= UpperAngle)
                    _limitState = 1; // At upper limit
                else
                    _limitState = 0; // Between limits
            }
            else
            {
                _limitImpulse = 0;
            }

            // Warm starting
            var P = _impulse;
            BodyA.Velocity -= P * mA;
            BodyA.AngularVelocity -= iA * (Vector2D.Cross(_rA, P) + _motorImpulse + _limitImpulse);

            if (BodyB != null)
            {
                BodyB.Velocity += P * mB;
                BodyB.AngularVelocity += iB * (Vector2D.Cross(_rB, P) + _motorImpulse + _limitImpulse);
            }
        }

        public override void SolveVelocity()
        {
            double iA = BodyA.InverseInertia;
            double iB = BodyB?.InverseInertia ?? 0;

            // Motor constraint
            if (EnableMotor && _limitState != 0)
            {
                double Cdot = (BodyB?.AngularVelocity ?? 0) - BodyA.AngularVelocity - MotorSpeed;
                double impulse = -_motorMass * Cdot;
                double oldImpulse = _motorImpulse;
                _motorImpulse = Math.Clamp(_motorImpulse + impulse, -MaxMotorTorque, MaxMotorTorque);
                impulse = _motorImpulse - oldImpulse;

                BodyA.AngularVelocity -= iA * impulse;
                if (BodyB != null)
                    BodyB.AngularVelocity += iB * impulse;
            }

            // Limit constraint
            if (EnableLimit && _limitState != 0)
            {
                double Cdot = (BodyB?.AngularVelocity ?? 0) - BodyA.AngularVelocity;
                double impulse = -_motorMass * Cdot;

                if (_limitState < 0)
                {
                    double old = _limitImpulse;
                    _limitImpulse = Math.Max(_limitImpulse + impulse, 0);
                    impulse = _limitImpulse - old;
                }
                else if (_limitState > 0)
                {
                    double old = _limitImpulse;
                    _limitImpulse = Math.Min(_limitImpulse + impulse, 0);
                    impulse = _limitImpulse - old;
                }

                BodyA.AngularVelocity -= iA * impulse;
                if (BodyB != null)
                    BodyB.AngularVelocity += iB * impulse;
            }

            // Point-to-point constraint
            var vA = BodyA.Velocity + Vector2D.Cross(BodyA.AngularVelocity, _rA);
            var vB = BodyB != null ? BodyB.Velocity + Vector2D.Cross(BodyB.AngularVelocity, _rB) : Vector2D.Zero;
            var Cdot2 = vB - vA;

            var impulse2 = SolveMatrix(-Cdot2);
            _impulse += impulse2;

            BodyA.Velocity -= impulse2 * BodyA.InverseMass;
            BodyA.AngularVelocity -= iA * Vector2D.Cross(_rA, impulse2);

            if (BodyB != null)
            {
                BodyB.Velocity += impulse2 * BodyB.InverseMass;
                BodyB.AngularVelocity += iB * Vector2D.Cross(_rB, impulse2);
            }
        }

        private Vector2D SolveMatrix(Vector2D b)
        {
            double det = _mass[0, 0] * _mass[1, 1] - _mass[0, 1] * _mass[1, 0];
            if (Math.Abs(det) < PhysicsConstants.Epsilon)
                return Vector2D.Zero;

            double invDet = 1.0 / det;
            return new Vector2D(
                invDet * (_mass[1, 1] * b.X - _mass[0, 1] * b.Y),
                invDet * (_mass[0, 0] * b.Y - _mass[1, 0] * b.X)
            );
        }

        public override bool SolvePosition()
        {
            var rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            var rB = BodyB != null ? Vector2D.Rotate(LocalAnchorB, BodyB.Rotation) : Vector2D.Zero;

            var C = (BodyB?.Position ?? Vector2D.Zero) + rB - BodyA.Position - rA;
            double error = C.Magnitude;

            double mA = BodyA.InverseMass;
            double mB = BodyB?.InverseMass ?? 0;
            double iA = BodyA.InverseInertia;
            double iB = BodyB?.InverseInertia ?? 0;

            var K = new double[2, 2];
            K[0, 0] = mA + mB + rA.Y * rA.Y * iA + rB.Y * rB.Y * iB;
            K[0, 1] = -rA.Y * rA.X * iA - rB.Y * rB.X * iB;
            K[1, 0] = K[0, 1];
            K[1, 1] = mA + mB + rA.X * rA.X * iA + rB.X * rB.X * iB;

            double det = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0];
            if (Math.Abs(det) < PhysicsConstants.Epsilon)
                return error < 0.005;

            double invDet = 1.0 / det;
            var impulse = new Vector2D(
                -invDet * (K[1, 1] * C.X - K[0, 1] * C.Y),
                -invDet * (K[0, 0] * C.Y - K[1, 0] * C.X)
            );

            BodyA.Position -= impulse * mA;
            BodyA.Rotation -= iA * Vector2D.Cross(rA, impulse);

            if (BodyB != null)
            {
                BodyB.Position += impulse * mB;
                BodyB.Rotation += iB * Vector2D.Cross(rB, impulse);
            }

            return error < 0.005;
        }

        public override Vector2D GetReactionForce(double invDt)
            => _impulse * invDt;

        public override double GetReactionTorque(double invDt)
            => (_motorImpulse + _limitImpulse) * invDt;
    }

    /// <summary>
    /// Prismatic joint - allows linear sliding along an axis.
    /// </summary>
    public class PrismaticJoint2D : Joint2D
    {
        /// <summary>Local anchor on body A.</summary>
        public Vector2D LocalAnchorA { get; set; }

        /// <summary>Local anchor on body B.</summary>
        public Vector2D LocalAnchorB { get; set; }

        /// <summary>Local axis on body A.</summary>
        public Vector2D LocalAxisA { get; set; }

        /// <summary>Reference angle.</summary>
        public double ReferenceAngle { get; set; }

        /// <summary>Whether translation limits are enabled.</summary>
        public bool EnableLimit { get; set; }

        /// <summary>Lower translation limit.</summary>
        public double LowerTranslation { get; set; }

        /// <summary>Upper translation limit.</summary>
        public double UpperTranslation { get; set; }

        /// <summary>Whether the motor is enabled.</summary>
        public bool EnableMotor { get; set; }

        /// <summary>Motor target speed.</summary>
        public double MotorSpeed { get; set; }

        /// <summary>Maximum motor force.</summary>
        public double MaxMotorForce { get; set; }

        // Solver temp
        private Vector2D _axis;
        private Vector2D _perp;
        private double _motorImpulse;
        private double _impulse;
        private double _limitState;

        public PrismaticJoint2D(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchor, Vector2D axis, string? id = null) : base(id)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = bodyA.WorldToLocal(anchor);
            LocalAnchorB = bodyB.WorldToLocal(anchor);
            LocalAxisA = bodyA.WorldToLocal(bodyA.Position + axis) - LocalAnchorA;
            LocalAxisA = LocalAxisA.Normalized;
            ReferenceAngle = bodyB.Rotation - bodyA.Rotation;
        }

        public override void PreSolve(double dt)
        {
            _axis = Vector2D.Rotate(LocalAxisA, BodyA.Rotation);
            _perp = _axis.Perpendicular;
        }

        public override void SolveVelocity()
        {
            if (BodyB == null) return;

            // Simplified velocity solve
            var vA = BodyA.Velocity;
            var vB = BodyB.Velocity;
            var vRel = vB - vA;

            // Motor
            if (EnableMotor)
            {
                double Cdot = Vector2D.Dot(vRel, _axis) - MotorSpeed;
                double impulse = -Cdot; // Simplified mass
                impulse = Math.Clamp(impulse, -MaxMotorForce, MaxMotorForce);
                _motorImpulse += impulse;

                var P = _axis * impulse;
                BodyA.Velocity -= P * BodyA.InverseMass;
                BodyB.Velocity += P * BodyB.InverseMass;
            }

            // Constraint along perpendicular
            double perpCdot = Vector2D.Dot(vRel, _perp);
            double perpImpulse = -perpCdot;
            _impulse += perpImpulse;

            var perpP = _perp * perpImpulse;
            BodyA.Velocity -= perpP * BodyA.InverseMass;
            BodyB.Velocity += perpP * BodyB.InverseMass;
        }

        public override bool SolvePosition()
        {
            // Simplified position solve
            return true;
        }

        public override Vector2D GetReactionForce(double invDt)
            => (_axis * _motorImpulse + _perp * _impulse) * invDt;

        public override double GetReactionTorque(double invDt) => 0;
    }

    /// <summary>
    /// Mouse joint - drags a body towards a target point.
    /// </summary>
    public class MouseJoint2D : Joint2D
    {
        /// <summary>Local anchor on the body.</summary>
        public Vector2D LocalAnchor { get; set; }

        /// <summary>World target point.</summary>
        public Vector2D Target { get; set; }

        /// <summary>Maximum force applied.</summary>
        public double MaxForce { get; set; } = 1000;

        /// <summary>Frequency for soft constraint.</summary>
        public double Frequency { get; set; } = 5;

        /// <summary>Damping ratio.</summary>
        public double DampingRatio { get; set; } = 0.7;

        // Solver temp
        private Vector2D _rA;
        private Vector2D _impulse;
        private double _gamma;
        private Vector2D _bias;
        private double[,] _mass = new double[2, 2];

        public MouseJoint2D(RigidBody2D body, Vector2D target, string? id = null) : base(id)
        {
            BodyA = body;
            Target = target;
            LocalAnchor = body.WorldToLocal(target);
        }

        public override void PreSolve(double dt)
        {
            _rA = Vector2D.Rotate(LocalAnchor, BodyA.Rotation);

            double mA = BodyA.InverseMass;
            double iA = BodyA.InverseInertia;

            // Soft constraint
            double omega = 2 * Math.PI * Frequency;
            double d = 2 * omega * DampingRatio;
            double k = omega * omega;

            double mass = mA + iA * _rA.MagnitudeSquared;
            mass = mass > 0 ? 1.0 / mass : 0;

            _gamma = dt * (d + dt * k);
            _gamma = _gamma > 0 ? 1.0 / _gamma : 0;

            var C = BodyA.Position + _rA - Target;
            _bias = C * (dt * k * _gamma);

            mass += _gamma;
            _mass[0, 0] = mass > 0 ? 1.0 / mass : 0;
            _mass[1, 1] = _mass[0, 0];
            _mass[0, 1] = 0;
            _mass[1, 0] = 0;

            // Warm starting
            _impulse *= 0.98;
            BodyA.Velocity += _impulse * mA;
            BodyA.AngularVelocity += iA * Vector2D.Cross(_rA, _impulse);
        }

        public override void SolveVelocity()
        {
            var vA = BodyA.Velocity + Vector2D.Cross(BodyA.AngularVelocity, _rA);
            var Cdot = vA + _bias + _impulse * _gamma;

            var impulse = new Vector2D(
                -_mass[0, 0] * Cdot.X,
                -_mass[1, 1] * Cdot.Y
            );

            var oldImpulse = _impulse;
            _impulse += impulse;

            double maxImpulse = MaxForce * (1.0 / 60.0); // Assume 60 Hz
            if (_impulse.MagnitudeSquared > maxImpulse * maxImpulse)
                _impulse = _impulse.Normalized * maxImpulse;

            impulse = _impulse - oldImpulse;

            BodyA.Velocity += impulse * BodyA.InverseMass;
            BodyA.AngularVelocity += BodyA.InverseInertia * Vector2D.Cross(_rA, impulse);
        }

        public override bool SolvePosition() => true;

        public override Vector2D GetReactionForce(double invDt)
            => _impulse * invDt;

        public override double GetReactionTorque(double invDt) => 0;
    }

    /// <summary>
    /// Weld joint - locks two bodies together.
    /// </summary>
    public class WeldJoint2D : Joint2D
    {
        /// <summary>Local anchor on body A.</summary>
        public Vector2D LocalAnchorA { get; set; }

        /// <summary>Local anchor on body B.</summary>
        public Vector2D LocalAnchorB { get; set; }

        /// <summary>Reference angle.</summary>
        public double ReferenceAngle { get; set; }

        /// <summary>Frequency for soft constraint.</summary>
        public double Frequency { get; set; }

        /// <summary>Damping ratio.</summary>
        public double DampingRatio { get; set; } = 0.7;

        // Solver temp
        private Vector2D _rA;
        private Vector2D _rB;
        private Vector2D _impulse;
        private double _angularImpulse;

        public WeldJoint2D(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchor, string? id = null) : base(id)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = bodyA.WorldToLocal(anchor);
            LocalAnchorB = bodyB.WorldToLocal(anchor);
            ReferenceAngle = bodyB.Rotation - bodyA.Rotation;
        }

        public override void PreSolve(double dt)
        {
            _rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            _rB = BodyB != null ? Vector2D.Rotate(LocalAnchorB, BodyB.Rotation) : Vector2D.Zero;
        }

        public override void SolveVelocity()
        {
            if (BodyB == null) return;

            // Angular constraint
            double wA = BodyA.AngularVelocity;
            double wB = BodyB.AngularVelocity;
            double Cdot = wB - wA;

            double iA = BodyA.InverseInertia;
            double iB = BodyB.InverseInertia;
            double angularMass = iA + iB;
            angularMass = angularMass > 0 ? 1.0 / angularMass : 0;

            double angularImpulse = -angularMass * Cdot;
            _angularImpulse += angularImpulse;

            BodyA.AngularVelocity -= iA * angularImpulse;
            BodyB.AngularVelocity += iB * angularImpulse;

            // Linear constraint
            var vA = BodyA.Velocity + Vector2D.Cross(BodyA.AngularVelocity, _rA);
            var vB = BodyB.Velocity + Vector2D.Cross(BodyB.AngularVelocity, _rB);
            var Cdot2 = vB - vA;

            double mA = BodyA.InverseMass;
            double mB = BodyB.InverseMass;
            double mass = mA + mB;
            mass = mass > 0 ? 1.0 / mass : 0;

            var impulse = -Cdot2 * mass;
            _impulse += impulse;

            BodyA.Velocity -= impulse * mA;
            BodyB.Velocity += impulse * mB;
        }

        public override bool SolvePosition()
        {
            if (BodyB == null) return true;

            var rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            var rB = Vector2D.Rotate(LocalAnchorB, BodyB.Rotation);

            var C = BodyB.Position + rB - BodyA.Position - rA;
            double angularC = BodyB.Rotation - BodyA.Rotation - ReferenceAngle;

            double error = C.Magnitude + Math.Abs(angularC);

            double mA = BodyA.InverseMass;
            double mB = BodyB.InverseMass;
            double iA = BodyA.InverseInertia;
            double iB = BodyB.InverseInertia;

            // Angular
            double angularMass = iA + iB;
            angularMass = angularMass > 0 ? 1.0 / angularMass : 0;
            double angularImpulse = -angularMass * angularC;

            BodyA.Rotation -= iA * angularImpulse;
            BodyB.Rotation += iB * angularImpulse;

            // Linear
            double mass = mA + mB;
            mass = mass > 0 ? 1.0 / mass : 0;
            var impulse = -C * mass;

            BodyA.Position -= impulse * mA;
            BodyB.Position += impulse * mB;

            return error < 0.005;
        }

        public override Vector2D GetReactionForce(double invDt)
            => _impulse * invDt;

        public override double GetReactionTorque(double invDt)
            => _angularImpulse * invDt;
    }

    /// <summary>
    /// Rope joint - enforces a maximum distance between two points.
    /// </summary>
    public class RopeJoint2D : Joint2D
    {
        /// <summary>Local anchor on body A.</summary>
        public Vector2D LocalAnchorA { get; set; }

        /// <summary>Local anchor on body B.</summary>
        public Vector2D LocalAnchorB { get; set; }

        /// <summary>Maximum length of the rope.</summary>
        public double MaxLength { get; set; }

        // Solver temp
        private Vector2D _u;
        private Vector2D _rA;
        private Vector2D _rB;
        private double _mass;
        private double _impulse;
        private double _length;

        public RopeJoint2D(RigidBody2D bodyA, RigidBody2D bodyB,
            Vector2D anchorA, Vector2D anchorB, double maxLength, string? id = null) : base(id)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = bodyA.WorldToLocal(anchorA);
            LocalAnchorB = bodyB.WorldToLocal(anchorB);
            MaxLength = maxLength;
        }

        public override void PreSolve(double dt)
        {
            if (BodyB == null) return;

            _rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            _rB = Vector2D.Rotate(LocalAnchorB, BodyB.Rotation);

            _u = BodyB.Position + _rB - BodyA.Position - _rA;
            _length = _u.Magnitude;

            if (_length > MaxLength)
            {
                _u /= _length;
            }
            else
            {
                _u = Vector2D.Zero;
                _impulse = 0;
                return;
            }

            double crA = Vector2D.Cross(_rA, _u);
            double crB = Vector2D.Cross(_rB, _u);
            double invMass = BodyA.InverseMass + BodyB.InverseMass +
                             BodyA.InverseInertia * crA * crA +
                             BodyB.InverseInertia * crB * crB;

            _mass = invMass > 0 ? 1.0 / invMass : 0;

            // Warm starting
            var P = _u * _impulse;
            BodyA.Velocity -= P * BodyA.InverseMass;
            BodyA.AngularVelocity -= BodyA.InverseInertia * Vector2D.Cross(_rA, P);
            BodyB.Velocity += P * BodyB.InverseMass;
            BodyB.AngularVelocity += BodyB.InverseInertia * Vector2D.Cross(_rB, P);
        }

        public override void SolveVelocity()
        {
            if (BodyB == null || _length <= MaxLength) return;

            var vA = BodyA.Velocity + Vector2D.Cross(BodyA.AngularVelocity, _rA);
            var vB = BodyB.Velocity + Vector2D.Cross(BodyB.AngularVelocity, _rB);

            double Cdot = Vector2D.Dot(_u, vB - vA);
            double impulse = -_mass * Cdot;

            // Only pull, never push
            double oldImpulse = _impulse;
            _impulse = Math.Min(_impulse + impulse, 0);
            impulse = _impulse - oldImpulse;

            var P = _u * impulse;
            BodyA.Velocity -= P * BodyA.InverseMass;
            BodyA.AngularVelocity -= BodyA.InverseInertia * Vector2D.Cross(_rA, P);
            BodyB.Velocity += P * BodyB.InverseMass;
            BodyB.AngularVelocity += BodyB.InverseInertia * Vector2D.Cross(_rB, P);
        }

        public override bool SolvePosition()
        {
            if (BodyB == null) return true;

            var rA = Vector2D.Rotate(LocalAnchorA, BodyA.Rotation);
            var rB = Vector2D.Rotate(LocalAnchorB, BodyB.Rotation);

            var u = BodyB.Position + rB - BodyA.Position - rA;
            double length = u.Magnitude;

            double C = length - MaxLength;
            if (C <= 0)
                return true;

            u /= length;

            double crA = Vector2D.Cross(rA, u);
            double crB = Vector2D.Cross(rB, u);
            double invMass = BodyA.InverseMass + BodyB.InverseMass +
                             BodyA.InverseInertia * crA * crA +
                             BodyB.InverseInertia * crB * crB;

            double impulse = invMass > 0 ? -C / invMass : 0;
            var P = u * impulse;

            BodyA.Position -= P * BodyA.InverseMass;
            BodyA.Rotation -= BodyA.InverseInertia * Vector2D.Cross(rA, P);
            BodyB.Position += P * BodyB.InverseMass;
            BodyB.Rotation += BodyB.InverseInertia * Vector2D.Cross(rB, P);

            return C < 0.005;
        }

        public override Vector2D GetReactionForce(double invDt)
            => _u * (_impulse * invDt);

        public override double GetReactionTorque(double invDt) => 0;
    }
}
