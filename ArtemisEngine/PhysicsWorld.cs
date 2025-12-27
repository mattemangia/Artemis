namespace ArtemisEngine;

public class PhysicsWorld
{
    public Vector2 Gravity { get; set; }
    public List<RigidBody> Bodies { get; private set; }

    private const int VelocityIterations = 6;
    private const int PositionIterations = 2;

    public PhysicsWorld(Vector2 gravity)
    {
        Gravity = gravity;
        Bodies = new List<RigidBody>();
    }

    public void AddBody(RigidBody body)
    {
        Bodies.Add(body);
    }

    public void RemoveBody(RigidBody body)
    {
        Bodies.Remove(body);
    }

    public void Step(float deltaTime)
    {
        // Apply gravity
        foreach (var body in Bodies)
        {
            if (!body.IsStatic && !body.IsKinematic)
            {
                body.ApplyForce(Gravity * body.Mass * deltaTime);
            }
        }

        // Update velocities and positions
        foreach (var body in Bodies)
        {
            body.Update(deltaTime);
        }

        // Collision detection and resolution
        for (int i = 0; i < Bodies.Count; i++)
        {
            for (int j = i + 1; j < Bodies.Count; j++)
            {
                var bodyA = Bodies[i];
                var bodyB = Bodies[j];

                if (bodyA.IsStatic && bodyB.IsStatic)
                    continue;

                if (DetectCollision(bodyA, bodyB, out var collision))
                {
                    ResolveCollision(bodyA, bodyB, collision);
                }
            }
        }
    }

    private bool DetectCollision(RigidBody bodyA, RigidBody bodyB, out Collision collision)
    {
        collision = new Collision();

        if (bodyA.Shape is CircleShape circleA && bodyB.Shape is CircleShape circleB)
        {
            return CircleVsCircle(bodyA, bodyB, circleA, circleB, out collision);
        }
        else if (bodyA.Shape is BoxShape boxA && bodyB.Shape is BoxShape boxB)
        {
            return BoxVsBox(bodyA, bodyB, boxA, boxB, out collision);
        }
        else if (bodyA.Shape is CircleShape circle && bodyB.Shape is BoxShape box)
        {
            return CircleVsBox(bodyA, bodyB, circle, box, out collision);
        }
        else if (bodyA.Shape is BoxShape box2 && bodyB.Shape is CircleShape circle2)
        {
            bool result = CircleVsBox(bodyB, bodyA, circle2, box2, out collision);
            if (result)
            {
                collision.Normal = -collision.Normal;
            }
            return result;
        }

        return false;
    }

    private bool CircleVsCircle(RigidBody bodyA, RigidBody bodyB, CircleShape circleA, CircleShape circleB, out Collision collision)
    {
        collision = new Collision();

        Vector2 delta = bodyB.Position - bodyA.Position;
        float distanceSquared = delta.LengthSquared;
        float radiusSum = circleA.Radius + circleB.Radius;

        if (distanceSquared >= radiusSum * radiusSum)
            return false;

        float distance = MathF.Sqrt(distanceSquared);

        collision.Normal = distance > 0 ? delta / distance : new Vector2(1, 0);
        collision.Penetration = radiusSum - distance;
        collision.ContactPoint = bodyA.Position + collision.Normal * circleA.Radius;

        return true;
    }

    private bool BoxVsBox(RigidBody bodyA, RigidBody bodyB, BoxShape boxA, BoxShape boxB, out Collision collision)
    {
        collision = new Collision();

        Vector2 delta = bodyB.Position - bodyA.Position;

        // Simple AABB collision for now (assumes no rotation or minimal rotation)
        float halfWidthA = boxA.Width / 2;
        float halfHeightA = boxA.Height / 2;
        float halfWidthB = boxB.Width / 2;
        float halfHeightB = boxB.Height / 2;

        float overlapX = (halfWidthA + halfWidthB) - MathF.Abs(delta.X);
        float overlapY = (halfHeightA + halfHeightB) - MathF.Abs(delta.Y);

        if (overlapX <= 0 || overlapY <= 0)
            return false;

        if (overlapX < overlapY)
        {
            collision.Normal = new Vector2(delta.X > 0 ? 1 : -1, 0);
            collision.Penetration = overlapX;
        }
        else
        {
            collision.Normal = new Vector2(0, delta.Y > 0 ? 1 : -1);
            collision.Penetration = overlapY;
        }

        collision.ContactPoint = bodyA.Position + delta * 0.5f;
        return true;
    }

    private bool CircleVsBox(RigidBody circleBody, RigidBody boxBody, CircleShape circle, BoxShape box, out Collision collision)
    {
        collision = new Collision();

        Vector2 delta = circleBody.Position - boxBody.Position;

        float halfW = box.Width / 2;
        float halfH = box.Height / 2;

        // Clamp circle center to box
        Vector2 closest = new Vector2(
            Math.Clamp(delta.X, -halfW, halfW),
            Math.Clamp(delta.Y, -halfH, halfH)
        );

        Vector2 localPoint = delta - closest;
        float distanceSquared = localPoint.LengthSquared;

        if (distanceSquared >= circle.Radius * circle.Radius)
            return false;

        float distance = MathF.Sqrt(distanceSquared);

        if (distance > 0)
        {
            collision.Normal = localPoint / distance;
        }
        else
        {
            // Circle center is inside the box
            if (MathF.Abs(delta.X) > MathF.Abs(delta.Y))
            {
                collision.Normal = new Vector2(delta.X > 0 ? 1 : -1, 0);
            }
            else
            {
                collision.Normal = new Vector2(0, delta.Y > 0 ? 1 : -1);
            }
        }

        collision.Penetration = circle.Radius - distance;
        collision.ContactPoint = circleBody.Position - collision.Normal * circle.Radius;

        return true;
    }

    private void ResolveCollision(RigidBody bodyA, RigidBody bodyB, Collision collision)
    {
        // Position correction
        float percent = 0.8f;
        float slop = 0.01f;
        Vector2 correction = collision.Normal *
            (Math.Max(collision.Penetration - slop, 0.0f) / (bodyA.InverseMass + bodyB.InverseMass)) * percent;

        if (!bodyA.IsStatic)
            bodyA.Position -= correction * bodyA.InverseMass;
        if (!bodyB.IsStatic)
            bodyB.Position += correction * bodyB.InverseMass;

        // Velocity resolution
        Vector2 relativeVelocity = bodyB.Velocity - bodyA.Velocity;
        float velocityAlongNormal = Vector2.Dot(relativeVelocity, collision.Normal);

        if (velocityAlongNormal > 0)
            return;

        float restitution = Math.Min(bodyA.Restitution, bodyB.Restitution);
        float impulseMagnitude = -(1 + restitution) * velocityAlongNormal;
        impulseMagnitude /= bodyA.InverseMass + bodyB.InverseMass;

        Vector2 impulse = collision.Normal * impulseMagnitude;

        if (!bodyA.IsStatic)
            bodyA.Velocity -= impulse * bodyA.InverseMass;
        if (!bodyB.IsStatic)
            bodyB.Velocity += impulse * bodyB.InverseMass;

        // Friction
        relativeVelocity = bodyB.Velocity - bodyA.Velocity;
        Vector2 tangent = relativeVelocity - collision.Normal * Vector2.Dot(relativeVelocity, collision.Normal);
        if (tangent.LengthSquared > 0.0001f)
        {
            tangent = tangent.Normalized;

            float frictionMagnitude = -Vector2.Dot(relativeVelocity, tangent);
            frictionMagnitude /= bodyA.InverseMass + bodyB.InverseMass;

            float mu = MathF.Sqrt(bodyA.Friction * bodyA.Friction + bodyB.Friction * bodyB.Friction);
            Vector2 frictionImpulse = tangent * Math.Clamp(frictionMagnitude, -impulseMagnitude * mu, impulseMagnitude * mu);

            if (!bodyA.IsStatic)
                bodyA.Velocity -= frictionImpulse * bodyA.InverseMass;
            if (!bodyB.IsStatic)
                bodyB.Velocity += frictionImpulse * bodyB.InverseMass;
        }
    }
}

public struct Collision
{
    public Vector2 Normal;
    public float Penetration;
    public Vector2 ContactPoint;
}
