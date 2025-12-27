using System;

namespace Artemis.Bodies
{
    /// <summary>
    /// Defines collision filtering using layers and masks.
    /// </summary>
    public struct CollisionFilter
    {
        /// <summary>
        /// The layer this body belongs to (bit flag).
        /// </summary>
        public uint Layer;

        /// <summary>
        /// The layers this body can collide with (bit mask).
        /// </summary>
        public uint Mask;

        /// <summary>
        /// The collision group (bodies in the same group don't collide).
        /// Set to 0 for no group.
        /// </summary>
        public int Group;

        /// <summary>
        /// Default filter that collides with everything.
        /// </summary>
        public static readonly CollisionFilter Default = new()
        {
            Layer = 0xFFFFFFFF,
            Mask = 0xFFFFFFFF,
            Group = 0
        };

        /// <summary>
        /// Filter that doesn't collide with anything.
        /// </summary>
        public static readonly CollisionFilter None = new()
        {
            Layer = 0,
            Mask = 0,
            Group = 0
        };

        /// <summary>
        /// Creates a collision filter for a specific layer.
        /// </summary>
        /// <param name="layer">Layer number (0-31).</param>
        /// <param name="collidesWithAll">Whether to collide with all layers.</param>
        public static CollisionFilter ForLayer(int layer, bool collidesWithAll = true)
        {
            uint layerBit = 1u << layer;
            return new CollisionFilter
            {
                Layer = layerBit,
                Mask = collidesWithAll ? 0xFFFFFFFF : layerBit,
                Group = 0
            };
        }

        /// <summary>
        /// Creates a collision filter that only collides with specific layers.
        /// </summary>
        /// <param name="layer">This body's layer.</param>
        /// <param name="collidesWithLayers">Layers to collide with.</param>
        public static CollisionFilter Create(int layer, params int[] collidesWithLayers)
        {
            uint mask = 0;
            foreach (int l in collidesWithLayers)
            {
                mask |= 1u << l;
            }

            return new CollisionFilter
            {
                Layer = 1u << layer,
                Mask = mask,
                Group = 0
            };
        }

        /// <summary>
        /// Checks if two filters should collide.
        /// </summary>
        public static bool ShouldCollide(CollisionFilter a, CollisionFilter b)
        {
            // Same non-zero group = no collision
            if (a.Group != 0 && a.Group == b.Group)
                return false;

            // Check layer/mask
            return (a.Layer & b.Mask) != 0 && (b.Layer & a.Mask) != 0;
        }

        /// <summary>
        /// Adds a layer to the mask.
        /// </summary>
        public CollisionFilter WithLayer(int layer)
        {
            Mask |= 1u << layer;
            return this;
        }

        /// <summary>
        /// Removes a layer from the mask.
        /// </summary>
        public CollisionFilter WithoutLayer(int layer)
        {
            Mask &= ~(1u << layer);
            return this;
        }

        /// <summary>
        /// Sets the collision group.
        /// </summary>
        public CollisionFilter InGroup(int group)
        {
            Group = group;
            return this;
        }
    }

    /// <summary>
    /// Predefined collision layers.
    /// </summary>
    public static class CollisionLayers
    {
        public const int Default = 0;
        public const int Static = 1;
        public const int Dynamic = 2;
        public const int Player = 3;
        public const int Enemy = 4;
        public const int Projectile = 5;
        public const int Trigger = 6;
        public const int Sensor = 7;
        public const int Water = 8;
        public const int Debris = 9;
        public const int Vehicle = 10;
        public const int Ragdoll = 11;

        // Reserved for user: 12-31
        public const int User1 = 12;
        public const int User2 = 13;
        public const int User3 = 14;
        public const int User4 = 15;
    }

    /// <summary>
    /// Defines how a body interacts with others.
    /// </summary>
    [Flags]
    public enum InteractionFlags
    {
        /// <summary>
        /// No special interaction.
        /// </summary>
        None = 0,

        /// <summary>
        /// Body is a trigger (detects overlaps but doesn't resolve collisions).
        /// </summary>
        IsTrigger = 1 << 0,

        /// <summary>
        /// Body is a sensor (like trigger but also affects forces like buoyancy zones).
        /// </summary>
        IsSensor = 1 << 1,

        /// <summary>
        /// Body generates continuous collision detection.
        /// </summary>
        UseCCD = 1 << 2,

        /// <summary>
        /// Body can push other bodies but cannot be pushed.
        /// </summary>
        Immovable = 1 << 3,

        /// <summary>
        /// Body ignores gravity.
        /// </summary>
        IgnoreGravity = 1 << 4,

        /// <summary>
        /// Body affects particles but particles don't affect body.
        /// </summary>
        AffectsParticlesOnly = 1 << 5,

        /// <summary>
        /// Body is affected by particles.
        /// </summary>
        AffectedByParticles = 1 << 6,

        /// <summary>
        /// Body generates collision events.
        /// </summary>
        GeneratesEvents = 1 << 7,

        /// <summary>
        /// Body is one-way (only blocks from one direction).
        /// </summary>
        OneWay = 1 << 8
    }
}
