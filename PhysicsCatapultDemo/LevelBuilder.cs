
namespace PhysicsCatapultDemo;

public class LevelBuilder
{
    public static List<GameObject> BuildLevel1()
    {
        var objects = new List<GameObject>();

        // Ground
        var groundShape = new BoxShape(50, 1);
        var groundBody = new RigidBody(new Vector2(25, -1), 0, groundShape, isStatic: true);
        objects.Add(new GameObject(groundBody, GameObjectType.Ground));

        // Build a structure with different materials
        float baseX = 35;
        float baseY = 1;

        // Bottom layer - Stone foundation
        objects.Add(CreateBlock(new Vector2(baseX - 2, baseY), 2, 2, Material.Stone));
        objects.Add(CreateBlock(new Vector2(baseX + 2, baseY), 2, 2, Material.Stone));

        // Second layer - Wood platform
        objects.Add(CreateBlock(new Vector2(baseX, baseY + 2.5f), 6, 0.5f, Material.Wood));

        // Third layer - Glass and wood boxes
        objects.Add(CreateBlock(new Vector2(baseX - 1.5f, baseY + 4), 1, 2, Material.Glass));
        objects.Add(CreateBlock(new Vector2(baseX + 1.5f, baseY + 4), 1, 2, Material.Glass));

        // Top layer - Wood beam
        objects.Add(CreateBlock(new Vector2(baseX, baseY + 6), 4, 0.5f, Material.Wood));

        // Add some extra blocks on the side
        objects.Add(CreateBlock(new Vector2(baseX + 5, baseY), 1.5f, 1.5f, Material.Wood));
        objects.Add(CreateBlock(new Vector2(baseX + 5, baseY + 2), 1.5f, 1.5f, Material.Glass));

        return objects;
    }

    public static List<GameObject> BuildLevel2()
    {
        var objects = new List<GameObject>();

        // Ground
        var groundShape = new BoxShape(50, 1);
        var groundBody = new RigidBody(new Vector2(25, -1), 0, groundShape, isStatic: true);
        objects.Add(new GameObject(groundBody, GameObjectType.Ground));

        // Build a tower
        float baseX = 35;
        float baseY = 1;

        for (int i = 0; i < 5; i++)
        {
            Material mat = i % 3 == 0 ? Material.Stone : (i % 3 == 1 ? Material.Wood : Material.Glass);
            objects.Add(CreateBlock(new Vector2(baseX, baseY + i * 2), 2, 1.5f, mat));
        }

        // Add side supports
        objects.Add(CreateBlock(new Vector2(baseX - 2.5f, baseY), 1, 3, Material.Stone));
        objects.Add(CreateBlock(new Vector2(baseX + 2.5f, baseY), 1, 3, Material.Stone));

        return objects;
    }

    public static List<GameObject> BuildLevel3()
    {
        var objects = new List<GameObject>();

        // Ground
        var groundShape = new BoxShape(50, 1);
        var groundBody = new RigidBody(new Vector2(25, -1), 0, groundShape, isStatic: true);
        objects.Add(new GameObject(groundBody, GameObjectType.Ground));

        // Build a pyramid
        float baseX = 35;
        float baseY = 1;

        // Bottom row - 5 blocks
        for (int i = 0; i < 5; i++)
        {
            objects.Add(CreateBlock(new Vector2(baseX - 4 + i * 2, baseY), 1.8f, 1.8f, Material.Stone));
        }

        // Second row - 4 blocks
        for (int i = 0; i < 4; i++)
        {
            objects.Add(CreateBlock(new Vector2(baseX - 3 + i * 2, baseY + 2.2f), 1.8f, 1.8f, Material.Wood));
        }

        // Third row - 3 blocks
        for (int i = 0; i < 3; i++)
        {
            objects.Add(CreateBlock(new Vector2(baseX - 2 + i * 2, baseY + 4.4f), 1.8f, 1.8f, Material.Glass));
        }

        // Fourth row - 2 blocks
        for (int i = 0; i < 2; i++)
        {
            objects.Add(CreateBlock(new Vector2(baseX - 1 + i * 2, baseY + 6.6f), 1.8f, 1.8f, Material.Wood));
        }

        // Top - 1 block
        objects.Add(CreateBlock(new Vector2(baseX, baseY + 8.8f), 1.8f, 1.8f, Material.Metal));

        return objects;
    }

    private static GameObject CreateBlock(Vector2 position, float width, float height, Material material)
    {
        var shape = new BoxShape(width * 0.5f, height * 0.5f);
        float volume = width * height;
        float mass = volume * material.Density;
        var body = new RigidBody(position, mass, shape);

        return new GameObject(body, GameObjectType.Block, material);
    }
}
