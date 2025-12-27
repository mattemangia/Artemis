namespace PhysicsCatapultDemo;

public class Material
{
    public string Name { get; set; }
    public float Density { get; set; }
    public float Strength { get; set; } // How much force needed to break
    public float Restitution { get; set; }
    public float Friction { get; set; }
    public ConsoleColor Color { get; set; }

    public Material(string name, float density, float strength, float restitution, float friction, ConsoleColor color)
    {
        Name = name;
        Density = density;
        Strength = strength;
        Restitution = restitution;
        Friction = friction;
        Color = color;
    }

    public static Material Wood => new("Wood", 0.6f, 50f, 0.3f, 0.4f, ConsoleColor.Yellow);
    public static Material Stone => new("Stone", 2.5f, 200f, 0.2f, 0.6f, ConsoleColor.Gray);
    public static Material Glass => new("Glass", 2.5f, 15f, 0.1f, 0.2f, ConsoleColor.Cyan);
    public static Material Metal => new("Metal", 7.8f, 500f, 0.4f, 0.3f, ConsoleColor.White);
}
