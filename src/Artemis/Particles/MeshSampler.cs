using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Numerics;

namespace Artemis.Particles
{
    /// <summary>
    /// A triangle mesh: vertex positions plus three vertex indices per triangle.
    /// </summary>
    public sealed class TriangleMesh
    {
        /// <summary>Vertex positions.</summary>
        public List<Vector3> Vertices { get; } = new();

        /// <summary>Triangle list: three indices into <see cref="Vertices"/> per triangle.</summary>
        public List<int> Indices { get; } = new();

        /// <summary>Number of triangles.</summary>
        public int TriangleCount => Indices.Count / 3;

        /// <summary>
        /// Axis-aligned bounds of the vertices.
        /// </summary>
        public (Vector3 Min, Vector3 Max) GetBounds()
        {
            if (Vertices.Count == 0)
                return (Vector3.Zero, Vector3.Zero);

            var min = Vertices[0];
            var max = Vertices[0];
            foreach (var v in Vertices)
            {
                min = Vector3.Min(min, v);
                max = Vector3.Max(max, v);
            }
            return (min, max);
        }

        /// <summary>
        /// Scales and moves the mesh in place: <c>v = v * scale + offset</c>.
        /// </summary>
        public TriangleMesh Transform(float scale, Vector3 offset)
        {
            for (int i = 0; i < Vertices.Count; i++)
                Vertices[i] = Vertices[i] * scale + offset;
            return this;
        }

        /// <summary>
        /// Loads the geometry of a Wavefront OBJ file. Only <c>v</c> and <c>f</c>
        /// records are read; polygons are split into triangle fans.
        /// </summary>
        public static TriangleMesh LoadObj(string path)
        {
            using var reader = new StreamReader(path);
            return ParseObj(reader);
        }

        /// <summary>
        /// Parses Wavefront OBJ text. See <see cref="LoadObj"/>.
        /// </summary>
        public static TriangleMesh ParseObj(string text)
        {
            using var reader = new StringReader(text);
            return ParseObj(reader);
        }

        /// <summary>
        /// Parses Wavefront OBJ text from a reader. See <see cref="LoadObj"/>.
        /// </summary>
        public static TriangleMesh ParseObj(TextReader reader)
        {
            var mesh = new TriangleMesh();
            var face = new List<int>(8);
            var separators = new[] { ' ', '\t' };

            string? line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.Length < 2)
                    continue;

                if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t'))
                {
                    var parts = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length >= 4)
                    {
                        mesh.Vertices.Add(new Vector3(
                            float.Parse(parts[1], CultureInfo.InvariantCulture),
                            float.Parse(parts[2], CultureInfo.InvariantCulture),
                            float.Parse(parts[3], CultureInfo.InvariantCulture)));
                    }
                }
                else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t'))
                {
                    var parts = line.Split(separators, StringSplitOptions.RemoveEmptyEntries);
                    face.Clear();
                    for (int i = 1; i < parts.Length; i++)
                    {
                        // "v", "v/vt", "v//vn" or "v/vt/vn"; indices are 1-based and
                        // negative ones count back from the last vertex read so far.
                        string token = parts[i];
                        int slash = token.IndexOf('/');
                        if (slash >= 0)
                            token = token.Substring(0, slash);

                        if (!int.TryParse(token, NumberStyles.Integer, CultureInfo.InvariantCulture, out int index) || index == 0)
                            continue;

                        index = index > 0 ? index - 1 : mesh.Vertices.Count + index;
                        if (index >= 0 && index < mesh.Vertices.Count)
                            face.Add(index);
                    }

                    for (int i = 2; i < face.Count; i++)
                    {
                        mesh.Indices.Add(face[0]);
                        mesh.Indices.Add(face[i - 1]);
                        mesh.Indices.Add(face[i]);
                    }
                }
            }

            return mesh;
        }

        /// <summary>
        /// Creates an axis-aligned box mesh (12 triangles).
        /// </summary>
        public static TriangleMesh CreateBox(Vector3 min, Vector3 max)
        {
            var mesh = new TriangleMesh();
            for (int i = 0; i < 8; i++)
            {
                mesh.Vertices.Add(new Vector3(
                    (i & 1) == 0 ? min.X : max.X,
                    (i & 2) == 0 ? min.Y : max.Y,
                    (i & 4) == 0 ? min.Z : max.Z));
            }

            mesh.Indices.AddRange(new[]
            {
                0, 2, 1, 1, 2, 3,  4, 5, 6, 5, 7, 6,
                0, 1, 4, 1, 5, 4,  2, 6, 3, 3, 6, 7,
                0, 4, 2, 2, 4, 6,  1, 3, 5, 3, 7, 5
            });
            return mesh;
        }

        /// <summary>
        /// Creates a closed surface of revolution around the Y axis from a profile of
        /// (radius, height) points, for example a vase or a drinking glass. With
        /// <paramref name="wallThickness"/> greater than zero the shape is hollow: the
        /// profile is the outer wall and an inner wall is added, leaving the top open
        /// and a solid base of the same thickness.
        /// </summary>
        public static TriangleMesh CreateLathe(IReadOnlyList<Vector2> profile, int segments = 32, float wallThickness = 0)
        {
            if (profile == null || profile.Count < 2)
                throw new ArgumentException("The profile needs at least two points.", nameof(profile));
            if (segments < 3)
                throw new ArgumentOutOfRangeException(nameof(segments));

            // Closed outline in the (radius, height) half plane, starting and ending
            // on the axis: outer wall bottom-up, then the inner wall top-down.
            var outline = new List<Vector2> { new(0, profile[0].Y) };
            outline.AddRange(profile);

            if (wallThickness > 0)
            {
                float floor = profile[0].Y + wallThickness;
                for (int i = profile.Count - 1; i >= 0; i--)
                {
                    float radius = Math.Max(profile[i].X - wallThickness, 0);
                    float height = Math.Max(profile[i].Y, floor);
                    if (radius > 0)
                        outline.Add(new Vector2(radius, height));
                }
                outline.Add(new Vector2(0, floor));
            }
            else
            {
                outline.Add(new Vector2(0, profile[profile.Count - 1].Y));
            }

            var mesh = new TriangleMesh();
            for (int p = 0; p < outline.Count; p++)
            {
                for (int s = 0; s < segments; s++)
                {
                    double angle = 2 * Math.PI * s / segments;
                    mesh.Vertices.Add(new Vector3(
                        outline[p].X * (float)Math.Cos(angle),
                        outline[p].Y,
                        outline[p].X * (float)Math.Sin(angle)));
                }
            }

            for (int p = 0; p + 1 < outline.Count; p++)
            {
                for (int s = 0; s < segments; s++)
                {
                    int a = p * segments + s;
                    int b = p * segments + (s + 1) % segments;
                    int c = a + segments;
                    int d = b + segments;

                    // Rings on the axis collapse to a point: skip the empty triangle.
                    if (outline[p].X > 0) { mesh.Indices.Add(a); mesh.Indices.Add(c); mesh.Indices.Add(b); }
                    if (outline[p + 1].X > 0) { mesh.Indices.Add(b); mesh.Indices.Add(c); mesh.Indices.Add(d); }
                }
            }

            return mesh;
        }
    }

    /// <summary>
    /// Turns triangle meshes into particle positions, so that a model can be rebuilt
    /// out of grains and then crumble, erode or shatter.
    /// </summary>
    public static class MeshSampler
    {
        /// <summary>
        /// Fills the volume enclosed by a closed mesh with points on a regular grid.
        /// <para>
        /// For every grid column a ray is cast through the mesh and the points between
        /// alternate surface crossings are kept (even-odd rule), so concave shapes and
        /// meshes with cavities are handled. The mesh should be watertight; open
        /// meshes give streaks, use <see cref="SampleSurface"/> for those.
        /// </para>
        /// </summary>
        /// <param name="mesh">The mesh to fill.</param>
        /// <param name="spacing">Distance between neighbouring points (one particle diameter).</param>
        /// <param name="maxPoints">Stop after this many points.</param>
        public static List<Vector3> FillVolume(TriangleMesh mesh, float spacing, int maxPoints = int.MaxValue)
        {
            if (mesh == null) throw new ArgumentNullException(nameof(mesh));
            if (!(spacing > 0)) throw new ArgumentOutOfRangeException(nameof(spacing));

            var points = new List<Vector3>();
            if (mesh.TriangleCount == 0)
                return points;

            var (min, max) = mesh.GetBounds();

            // Sample rays run along X through the centres of the (y, z) cells. A tiny
            // irrational-looking offset keeps them from passing exactly through the
            // vertices and edges of axis-aligned models, where a crossing could be
            // counted twice or not at all.
            float originY = min.Y + spacing * 0.5f + spacing * 0.00037f;
            float originZ = min.Z + spacing * 0.5f + spacing * 0.00061f;
            float originX = min.X + spacing * 0.5f;

            int ny = Math.Max(1, (int)Math.Ceiling((max.Y - min.Y) / spacing));
            int nz = Math.Max(1, (int)Math.Ceiling((max.Z - min.Z) / spacing));
            int nx = Math.Max(1, (int)Math.Ceiling((max.X - min.X) / spacing));

            // Triangles bucketed by the columns their (y, z) bounds cover.
            var columns = new List<int>?[ny * nz];
            var vertices = mesh.Vertices;
            var indices = mesh.Indices;

            for (int t = 0; t < indices.Count; t += 3)
            {
                var a = vertices[indices[t]];
                var b = vertices[indices[t + 1]];
                var c = vertices[indices[t + 2]];

                int y0 = Math.Max(0, (int)Math.Floor((Math.Min(a.Y, Math.Min(b.Y, c.Y)) - originY) / spacing));
                int y1 = Math.Min(ny - 1, (int)Math.Ceiling((Math.Max(a.Y, Math.Max(b.Y, c.Y)) - originY) / spacing));
                int z0 = Math.Max(0, (int)Math.Floor((Math.Min(a.Z, Math.Min(b.Z, c.Z)) - originZ) / spacing));
                int z1 = Math.Min(nz - 1, (int)Math.Ceiling((Math.Max(a.Z, Math.Max(b.Z, c.Z)) - originZ) / spacing));

                for (int z = z0; z <= z1; z++)
                {
                    for (int y = y0; y <= y1; y++)
                        (columns[z * ny + y] ??= new List<int>(4)).Add(t);
                }
            }

            var crossings = new List<float>(16);

            for (int z = 0; z < nz; z++)
            {
                float rayZ = originZ + z * spacing;
                for (int y = 0; y < ny; y++)
                {
                    var bucket = columns[z * ny + y];
                    if (bucket == null)
                        continue;

                    float rayY = originY + y * spacing;
                    crossings.Clear();

                    foreach (int t in bucket)
                    {
                        if (RayCrossing(vertices[indices[t]], vertices[indices[t + 1]], vertices[indices[t + 2]],
                                rayY, rayZ, out float x))
                            crossings.Add(x);
                    }

                    if (crossings.Count < 2)
                        continue;

                    crossings.Sort();

                    for (int i = 0; i + 1 < crossings.Count; i += 2)
                    {
                        int from = Math.Max(0, (int)Math.Ceiling((crossings[i] - originX) / spacing));
                        int to = Math.Min(nx - 1, (int)Math.Floor((crossings[i + 1] - originX) / spacing));
                        for (int x = from; x <= to; x++)
                        {
                            points.Add(new Vector3(originX + x * spacing, rayY, rayZ));
                            if (points.Count >= maxPoints)
                                return points;
                        }
                    }
                }
            }

            return points;
        }

        /// <summary>
        /// Where the line through (y, z) parallel to the X axis crosses a triangle.
        /// </summary>
        private static bool RayCrossing(Vector3 a, Vector3 b, Vector3 c, float y, float z, out float x)
        {
            x = 0;

            // Barycentric coordinates of (y, z) in the triangle projected on the YZ plane.
            double d = (double)(b.Y - a.Y) * (c.Z - a.Z) - (double)(c.Y - a.Y) * (b.Z - a.Z);
            if (Math.Abs(d) < 1e-20)
                return false; // triangle is edge-on to the ray

            double u = ((double)(y - a.Y) * (c.Z - a.Z) - (double)(c.Y - a.Y) * (z - a.Z)) / d;
            double v = ((double)(b.Y - a.Y) * (z - a.Z) - (double)(y - a.Y) * (b.Z - a.Z)) / d;
            if (u < 0 || v < 0 || u + v > 1)
                return false;

            x = (float)(a.X + u * (b.X - a.X) + v * (c.X - a.X));
            return true;
        }

        /// <summary>
        /// Covers the surface of a mesh with points roughly <paramref name="spacing"/>
        /// apart. Works for open meshes as well; use it for shells, or when a model is
        /// too large to fill solid.
        /// </summary>
        public static List<Vector3> SampleSurface(TriangleMesh mesh, float spacing, int maxPoints = int.MaxValue)
        {
            if (mesh == null) throw new ArgumentNullException(nameof(mesh));
            if (!(spacing > 0)) throw new ArgumentOutOfRangeException(nameof(spacing));

            var points = new List<Vector3>();
            var occupied = new HashSet<(int, int, int)>();
            float inv = 1f / spacing;
            var vertices = mesh.Vertices;
            var indices = mesh.Indices;

            for (int t = 0; t < indices.Count; t += 3)
            {
                var a = vertices[indices[t]];
                var b = vertices[indices[t + 1]];
                var c = vertices[indices[t + 2]];

                float longest = Math.Max((b - a).Length(), Math.Max((c - a).Length(), (c - b).Length()));
                int steps = Math.Max(1, (int)Math.Ceiling(longest / (spacing * 0.5f)));

                for (int i = 0; i <= steps; i++)
                {
                    for (int j = 0; i + j <= steps; j++)
                    {
                        float u = (float)i / steps, v = (float)j / steps;
                        var p = a + (b - a) * u + (c - a) * v;

                        // One point per grid cell keeps the density even where
                        // triangles are small or overlap.
                        var cell = ((int)MathF.Floor(p.X * inv), (int)MathF.Floor(p.Y * inv), (int)MathF.Floor(p.Z * inv));
                        if (!occupied.Add(cell))
                            continue;

                        points.Add(new Vector3((cell.Item1 + 0.5f) * spacing, (cell.Item2 + 0.5f) * spacing, (cell.Item3 + 0.5f) * spacing));
                        if (points.Count >= maxPoints)
                            return points;
                    }
                }
            }

            return points;
        }
    }
}
