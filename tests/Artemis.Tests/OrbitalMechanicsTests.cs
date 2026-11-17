using Artemis.Core;
using Artemis.Forces;

namespace Artemis.Tests;

public class OrbitalMechanicsTests
{
    private static (OrbitalMechanics System, CelestialBody Earth, CelestialBody Satellite) Setup(bool nBody)
    {
        var system = new OrbitalMechanics { UseNBody = nBody };
        var earth = OrbitalMechanics.CreateEarth();
        earth.IsFixed = true;
        earth.Position = Vector3D.Zero;
        earth.Velocity = Vector3D.Zero;
        var satellite = new CelestialBody { Mass = 1000, Radius = 1 };
        system.AddBody(earth);
        system.AddBody(satellite);
        return (system, earth, satellite);
    }

    [Fact]
    public void N_body_circular_orbit_keeps_its_radius_and_has_the_kepler_period()
    {
        var (system, earth, satellite) = Setup(nBody: true);
        const double r = 7.0e6;
        system.SetCircularOrbit(satellite, earth, r);
        double expectedPeriod = 2 * Math.PI * Math.Sqrt(r * r * r / earth.Mu);

        double minR = double.MaxValue, maxR = 0, period = 0, t = 0;
        double previousY = satellite.Position.Y;
        while (t < 1.5 * expectedPeriod)
        {
            system.Update(1);
            t += 1;
            double radius = satellite.Position.Magnitude;
            minR = Math.Min(minR, radius);
            maxR = Math.Max(maxR, radius);
            if (period == 0 && t > 0.5 * expectedPeriod && previousY < 0 && satellite.Position.Y >= 0)
                period = t;
            previousY = satellite.Position.Y;
        }

        Assert.InRange(minR / r, 0.9999, 1.0001);
        Assert.InRange(maxR / r, 0.9999, 1.0001);
        Assert.InRange(period / expectedPeriod, 0.999, 1.001);
    }

    [Fact]
    public void Kepler_orbit_set_after_the_simulation_started_begins_at_periapsis_and_closes()
    {
        var (system, earth, satellite) = Setup(nBody: false);
        system.SetCircularOrbit(satellite, earth, 7.0e6);
        for (int i = 0; i < 1000; i++) system.Update(1); // simulation time is no longer zero

        system.SetEllipticalOrbit(satellite, earth, 7.0e6, 14.0e6);
        var elements = system.CalculateOrbitalElements(satellite, earth);
        Assert.InRange(elements.SemiMajorAxis, 1.0499e7, 1.0501e7);
        Assert.Equal(1.0 / 3.0, elements.Eccentricity, 6);

        system.Update(1);
        Assert.InRange(satellite.Position.Magnitude, 7.0e6, 7.001e6); // still next to periapsis

        double minR = double.MaxValue, maxR = 0;
        int steps = (int)Math.Round(elements.Period) - 1;
        for (int i = 0; i < steps; i++)
        {
            system.Update(1);
            double radius = satellite.Position.Magnitude;
            minR = Math.Min(minR, radius);
            maxR = Math.Max(maxR, radius);
        }

        Assert.InRange(minR, 6.999e6, 7.001e6);
        Assert.InRange(maxR, 13.999e6, 14.001e6);
        Assert.InRange(satellite.Position.Magnitude, 7.0e6, 7.001e6); // back at periapsis after one period
    }
}
