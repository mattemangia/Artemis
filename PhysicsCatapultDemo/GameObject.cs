using ArtemisEngine;

namespace PhysicsCatapultDemo;

public class GameObject
{
    public RigidBody Body { get; set; }
    public Material? Material { get; set; }
    public float Health { get; set; }
    public float MaxHealth { get; set; }
    public bool IsDestroyed { get; set; }
    public GameObjectType Type { get; set; }

    public GameObject(RigidBody body, GameObjectType type, Material? material = null)
    {
        Body = body;
        Type = type;
        Material = material;

        if (material != null)
        {
            MaxHealth = material.Strength;
            Health = MaxHealth;
            body.Restitution = material.Restitution;
            body.Friction = material.Friction;
        }
        else
        {
            MaxHealth = 100;
            Health = MaxHealth;
        }
    }

    public void TakeDamage(float damage)
    {
        Health -= damage;
        if (Health <= 0)
        {
            IsDestroyed = true;
        }
    }
}

public enum GameObjectType
{
    Projectile,
    Block,
    Ground,
    Catapult
}
