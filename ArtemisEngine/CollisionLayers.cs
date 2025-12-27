namespace ArtemisEngine;

/// <summary>
/// Collision layers system for filtering what collides with what
/// Uses bit masking for efficient collision filtering
/// </summary>
public static class CollisionLayers
{
    public const int Default = 1 << 0;      // 0001
    public const int Static = 1 << 1;       // 0010
    public const int Projectile = 1 << 2;   // 0100
    public const int Enemy = 1 << 3;        // 1000
    public const int Player = 1 << 4;       // 0001 0000
    public const int Trigger = 1 << 5;      // 0010 0000
    public const int Ground = 1 << 6;       // 0100 0000
    public const int Debris = 1 << 7;       // 1000 0000

    public const int Everything = ~0; // All bits set

    /// <summary>
    /// Check if two layers should collide
    /// </summary>
    public static bool ShouldCollide(int layerA, int maskA, int layerB, int maskB)
    {
        return (layerA & maskB) != 0 && (layerB & maskA) != 0;
    }
}

public class CollisionFilter
{
    public int Layer { get; set; } = CollisionLayers.Default;
    public int Mask { get; set; } = CollisionLayers.Everything;

    public CollisionFilter() { }

    public CollisionFilter(int layer, int mask)
    {
        Layer = layer;
        Mask = mask;
    }

    public bool CanCollideWith(CollisionFilter other)
    {
        return CollisionLayers.ShouldCollide(Layer, Mask, other.Layer, other.Mask);
    }
}
