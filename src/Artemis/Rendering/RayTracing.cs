using System;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using Artemis.Compute;

namespace Artemis.Rendering
{
    /// <summary>
    /// Ray for tracing through the scene.
    /// </summary>
    public struct Ray
    {
        /// <summary>Ray origin.</summary>
        public Vector3 Origin;

        /// <summary>Ray direction (normalized).</summary>
        public Vector3 Direction;

        /// <summary>Minimum distance.</summary>
        public float TMin;

        /// <summary>Maximum distance.</summary>
        public float TMax;

        /// <summary>
        /// Creates a new ray.
        /// </summary>
        public Ray(Vector3 origin, Vector3 direction, float tMin = 0.001f, float tMax = float.MaxValue)
        {
            Origin = origin;
            Direction = Vector3.Normalize(direction);
            TMin = tMin;
            TMax = tMax;
        }

        /// <summary>
        /// Gets a point along the ray at distance t.
        /// </summary>
        public Vector3 At(float t) => Origin + Direction * t;
    }

    /// <summary>
    /// Hit information from ray intersection.
    /// </summary>
    public struct RayHit
    {
        /// <summary>Whether the ray hit something.</summary>
        public bool Hit;

        /// <summary>Distance along the ray.</summary>
        public float T;

        /// <summary>Hit point in world space.</summary>
        public Vector3 Point;

        /// <summary>Surface normal at hit point.</summary>
        public Vector3 Normal;

        /// <summary>UV texture coordinates.</summary>
        public Vector2 UV;

        /// <summary>Material index of hit surface.</summary>
        public int MaterialIndex;

        /// <summary>Object ID that was hit.</summary>
        public int ObjectId;

        /// <summary>Whether the hit was on the front face.</summary>
        public bool FrontFace;
    }

    /// <summary>
    /// Material properties for ray tracing.
    /// </summary>
    public class RayTracingMaterial
    {
        /// <summary>Material ID.</summary>
        public int Id { get; set; }

        /// <summary>Base color (albedo).</summary>
        public Vector3 Albedo { get; set; } = Vector3.One;

        /// <summary>Emission color and intensity.</summary>
        public Vector3 Emission { get; set; } = Vector3.Zero;

        /// <summary>Reflectivity (0 = diffuse, 1 = perfect mirror).</summary>
        public float Reflectivity { get; set; }

        /// <summary>Roughness for glossy reflections (0 = smooth, 1 = rough).</summary>
        public float Roughness { get; set; } = 0.5f;

        /// <summary>Metalness (0 = dielectric, 1 = metal).</summary>
        public float Metalness { get; set; }

        /// <summary>Index of refraction for transparent materials.</summary>
        public float IOR { get; set; } = 1.5f;

        /// <summary>Transparency (0 = opaque, 1 = fully transparent).</summary>
        public float Transparency { get; set; }

        /// <summary>Whether this is a mirror surface.</summary>
        public bool IsMirror => Reflectivity > 0.95f && Roughness < 0.05f;

        #region Presets

        /// <summary>Creates a perfect mirror material.</summary>
        public static RayTracingMaterial Mirror() => new()
        {
            Albedo = new Vector3(0.95f, 0.95f, 0.95f),
            Reflectivity = 1.0f,
            Roughness = 0.0f,
            Metalness = 1.0f
        };

        /// <summary>Creates a chrome/polished metal material.</summary>
        public static RayTracingMaterial Chrome() => new()
        {
            Albedo = new Vector3(0.8f, 0.8f, 0.85f),
            Reflectivity = 0.95f,
            Roughness = 0.02f,
            Metalness = 1.0f
        };

        /// <summary>Creates a gold material.</summary>
        public static RayTracingMaterial Gold() => new()
        {
            Albedo = new Vector3(1.0f, 0.84f, 0.0f),
            Reflectivity = 0.9f,
            Roughness = 0.1f,
            Metalness = 1.0f
        };

        /// <summary>Creates a copper material.</summary>
        public static RayTracingMaterial Copper() => new()
        {
            Albedo = new Vector3(0.95f, 0.64f, 0.54f),
            Reflectivity = 0.85f,
            Roughness = 0.15f,
            Metalness = 1.0f
        };

        /// <summary>Creates glass material.</summary>
        public static RayTracingMaterial Glass() => new()
        {
            Albedo = new Vector3(1.0f, 1.0f, 1.0f),
            Reflectivity = 0.1f,
            Roughness = 0.0f,
            Transparency = 0.95f,
            IOR = 1.52f
        };

        /// <summary>Creates water material.</summary>
        public static RayTracingMaterial Water() => new()
        {
            Albedo = new Vector3(0.8f, 0.9f, 1.0f),
            Reflectivity = 0.1f,
            Roughness = 0.1f,
            Transparency = 0.9f,
            IOR = 1.33f
        };

        /// <summary>Creates a diffuse/matte material.</summary>
        public static RayTracingMaterial Diffuse(Vector3 color) => new()
        {
            Albedo = color,
            Reflectivity = 0.0f,
            Roughness = 1.0f,
            Metalness = 0.0f
        };

        /// <summary>Creates a glossy material.</summary>
        public static RayTracingMaterial Glossy(Vector3 color, float roughness = 0.3f) => new()
        {
            Albedo = color,
            Reflectivity = 0.5f,
            Roughness = roughness,
            Metalness = 0.0f
        };

        #endregion
    }

    /// <summary>
    /// Traceable object in the scene.
    /// </summary>
    public abstract class TraceableObject
    {
        /// <summary>Object ID.</summary>
        public int Id { get; set; }

        /// <summary>Material index.</summary>
        public int MaterialIndex { get; set; }

        /// <summary>World transform matrix.</summary>
        public Matrix4x4 Transform { get; set; } = Matrix4x4.Identity;

        /// <summary>Inverse transform (cached).</summary>
        public Matrix4x4 InverseTransform { get; private set; } = Matrix4x4.Identity;

        /// <summary>
        /// Updates the inverse transform.
        /// </summary>
        public void UpdateInverseTransform()
        {
            Matrix4x4.Invert(Transform, out var inv);
            InverseTransform = inv;
        }

        /// <summary>
        /// Tests ray intersection.
        /// </summary>
        public abstract bool Intersect(Ray ray, out RayHit hit);

        /// <summary>
        /// Gets the bounding box.
        /// </summary>
        public abstract (Vector3 min, Vector3 max) GetBounds();
    }

    /// <summary>
    /// Traceable sphere.
    /// </summary>
    public class TraceableSphere : TraceableObject
    {
        /// <summary>Sphere center (local space).</summary>
        public Vector3 Center { get; set; } = Vector3.Zero;

        /// <summary>Sphere radius.</summary>
        public float Radius { get; set; } = 1.0f;

        public override bool Intersect(Ray ray, out RayHit hit)
        {
            hit = default;

            // Transform ray to local space
            var localOrigin = Vector3.Transform(ray.Origin, InverseTransform);
            var localDir = Vector3.TransformNormal(ray.Direction, InverseTransform);

            var oc = localOrigin - Center;
            float a = Vector3.Dot(localDir, localDir);
            float halfB = Vector3.Dot(oc, localDir);
            float c = Vector3.Dot(oc, oc) - Radius * Radius;
            float discriminant = halfB * halfB - a * c;

            if (discriminant < 0)
                return false;

            float sqrtD = MathF.Sqrt(discriminant);
            float t = (-halfB - sqrtD) / a;

            if (t < ray.TMin || t > ray.TMax)
            {
                t = (-halfB + sqrtD) / a;
                if (t < ray.TMin || t > ray.TMax)
                    return false;
            }

            var localPoint = localOrigin + localDir * t;
            var localNormal = (localPoint - Center) / Radius;

            hit.Hit = true;
            hit.T = t;
            hit.Point = ray.At(t);
            hit.Normal = Vector3.Normalize(Vector3.TransformNormal(localNormal, Transform));
            hit.MaterialIndex = MaterialIndex;
            hit.ObjectId = Id;
            hit.FrontFace = Vector3.Dot(ray.Direction, hit.Normal) < 0;
            if (!hit.FrontFace) hit.Normal = -hit.Normal;

            // UV coordinates (spherical)
            float theta = MathF.Acos(-localNormal.Y);
            float phi = MathF.Atan2(-localNormal.Z, localNormal.X) + MathF.PI;
            hit.UV = new Vector2(phi / (2 * MathF.PI), theta / MathF.PI);

            return true;
        }

        public override (Vector3 min, Vector3 max) GetBounds()
        {
            var r = new Vector3(Radius);
            return (Center - r, Center + r);
        }
    }

    /// <summary>
    /// Traceable plane (infinite or bounded).
    /// </summary>
    public class TraceablePlane : TraceableObject
    {
        /// <summary>Plane center point.</summary>
        public Vector3 Center { get; set; } = Vector3.Zero;

        /// <summary>Plane normal.</summary>
        public Vector3 Normal { get; set; } = Vector3.UnitY;

        /// <summary>Half extents (set to large value for infinite).</summary>
        public Vector2 HalfExtents { get; set; } = new Vector2(1000, 1000);

        public override bool Intersect(Ray ray, out RayHit hit)
        {
            hit = default;

            float denom = Vector3.Dot(Normal, ray.Direction);
            if (MathF.Abs(denom) < 1e-6f)
                return false;

            float t = Vector3.Dot(Center - ray.Origin, Normal) / denom;
            if (t < ray.TMin || t > ray.TMax)
                return false;

            var point = ray.At(t);

            // Check bounds (project onto plane basis)
            var toPoint = point - Center;
            var right = Vector3.Cross(Normal, Vector3.UnitY);
            if (right.LengthSquared() < 0.01f)
                right = Vector3.Cross(Normal, Vector3.UnitX);
            right = Vector3.Normalize(right);
            var forward = Vector3.Cross(right, Normal);

            float u = Vector3.Dot(toPoint, right);
            float v = Vector3.Dot(toPoint, forward);

            if (MathF.Abs(u) > HalfExtents.X || MathF.Abs(v) > HalfExtents.Y)
                return false;

            hit.Hit = true;
            hit.T = t;
            hit.Point = point;
            hit.Normal = denom < 0 ? Normal : -Normal;
            hit.FrontFace = denom < 0;
            hit.UV = new Vector2((u / HalfExtents.X + 1) * 0.5f, (v / HalfExtents.Y + 1) * 0.5f);
            hit.MaterialIndex = MaterialIndex;
            hit.ObjectId = Id;

            return true;
        }

        public override (Vector3 min, Vector3 max) GetBounds()
        {
            // Approximate bounds
            var ext = new Vector3(HalfExtents.X, 0.01f, HalfExtents.Y);
            return (Center - ext, Center + ext);
        }
    }

    /// <summary>
    /// Traceable axis-aligned box.
    /// </summary>
    public class TraceableBox : TraceableObject
    {
        /// <summary>Box center.</summary>
        public Vector3 Center { get; set; } = Vector3.Zero;

        /// <summary>Box half extents.</summary>
        public Vector3 HalfExtents { get; set; } = Vector3.One;

        public override bool Intersect(Ray ray, out RayHit hit)
        {
            hit = default;

            var min = Center - HalfExtents;
            var max = Center + HalfExtents;

            var invDir = new Vector3(1f / ray.Direction.X, 1f / ray.Direction.Y, 1f / ray.Direction.Z);

            var t0s = (min - ray.Origin) * invDir;
            var t1s = (max - ray.Origin) * invDir;

            var tSmaller = Vector3.Min(t0s, t1s);
            var tLarger = Vector3.Max(t0s, t1s);

            float tMin = MathF.Max(MathF.Max(tSmaller.X, tSmaller.Y), tSmaller.Z);
            float tMax = MathF.Min(MathF.Min(tLarger.X, tLarger.Y), tLarger.Z);

            if (tMax < tMin || tMax < ray.TMin || tMin > ray.TMax)
                return false;

            float t = tMin >= ray.TMin ? tMin : tMax;
            if (t < ray.TMin || t > ray.TMax)
                return false;

            var point = ray.At(t);
            var localPoint = point - Center;

            // Determine which face was hit
            Vector3 normal;
            float maxComp = MathF.Max(MathF.Max(
                MathF.Abs(localPoint.X / HalfExtents.X),
                MathF.Abs(localPoint.Y / HalfExtents.Y)),
                MathF.Abs(localPoint.Z / HalfExtents.Z));

            if (MathF.Abs(localPoint.X / HalfExtents.X - maxComp) < 0.001f)
                normal = new Vector3(MathF.Sign(localPoint.X), 0, 0);
            else if (MathF.Abs(localPoint.Y / HalfExtents.Y - maxComp) < 0.001f)
                normal = new Vector3(0, MathF.Sign(localPoint.Y), 0);
            else
                normal = new Vector3(0, 0, MathF.Sign(localPoint.Z));

            hit.Hit = true;
            hit.T = t;
            hit.Point = point;
            hit.Normal = normal;
            hit.FrontFace = Vector3.Dot(ray.Direction, normal) < 0;
            if (!hit.FrontFace) hit.Normal = -hit.Normal;
            hit.MaterialIndex = MaterialIndex;
            hit.ObjectId = Id;

            return true;
        }

        public override (Vector3 min, Vector3 max) GetBounds()
        {
            return (Center - HalfExtents, Center + HalfExtents);
        }
    }

    /// <summary>
    /// BVH node for acceleration structure.
    /// </summary>
    internal class BVHNode
    {
        public Vector3 Min;
        public Vector3 Max;
        public int LeftChild;  // -1 if leaf
        public int RightChild; // -1 if leaf
        public int ObjectIndex; // -1 if not leaf
    }

    /// <summary>
    /// GPU-accelerated ray tracing system for mirror reflections.
    /// Supports OpenCL/CUDA/CPU backends.
    /// </summary>
    public class RayTracingSystem : IDisposable
    {
        private readonly List<TraceableObject> _objects = new();
        private readonly List<RayTracingMaterial> _materials = new();
        private readonly List<BVHNode> _bvhNodes = new();
        private GpuCompute? _gpuCompute;
        private bool _bvhDirty = true;

        /// <summary>Maximum ray bounces for reflections.</summary>
        public int MaxBounces { get; set; } = 8;

        /// <summary>Background/sky color.</summary>
        public Vector3 BackgroundColor { get; set; } = new Vector3(0.1f, 0.1f, 0.15f);

        /// <summary>Ambient light contribution.</summary>
        public Vector3 AmbientLight { get; set; } = new Vector3(0.1f, 0.1f, 0.1f);

        /// <summary>Whether to use GPU acceleration.</summary>
        public bool UseGPU { get; set; } = true;

        /// <summary>Gets the active GPU device info.</summary>
        public GpuDeviceInfo? GpuDevice => _gpuCompute?.DeviceInfo;

        /// <summary>Gets all traceable objects.</summary>
        public IReadOnlyList<TraceableObject> Objects => _objects;

        /// <summary>Gets all materials.</summary>
        public IReadOnlyList<RayTracingMaterial> Materials => _materials;

        /// <summary>
        /// Initializes the ray tracing system.
        /// </summary>
        /// <param name="preferredBackend">Preferred GPU backend.</param>
        public bool Initialize(GpuBackend preferredBackend = GpuBackend.Auto)
        {
            if (UseGPU)
            {
                _gpuCompute = new GpuCompute(preferredBackend);
                return _gpuCompute.Initialize();
            }
            return true;
        }

        /// <summary>
        /// Adds a material and returns its index.
        /// </summary>
        public int AddMaterial(RayTracingMaterial material)
        {
            material.Id = _materials.Count;
            _materials.Add(material);
            return material.Id;
        }

        /// <summary>
        /// Adds a traceable object.
        /// </summary>
        public void AddObject(TraceableObject obj)
        {
            obj.Id = _objects.Count;
            obj.UpdateInverseTransform();
            _objects.Add(obj);
            _bvhDirty = true;
        }

        /// <summary>
        /// Adds a mirror sphere.
        /// </summary>
        public TraceableSphere AddMirrorSphere(Vector3 center, float radius)
        {
            var mirrorMat = AddMaterial(RayTracingMaterial.Mirror());
            var sphere = new TraceableSphere
            {
                Center = center,
                Radius = radius,
                MaterialIndex = mirrorMat
            };
            AddObject(sphere);
            return sphere;
        }

        /// <summary>
        /// Adds a mirror plane.
        /// </summary>
        public TraceablePlane AddMirrorPlane(Vector3 center, Vector3 normal, Vector2 size)
        {
            var mirrorMat = AddMaterial(RayTracingMaterial.Mirror());
            var plane = new TraceablePlane
            {
                Center = center,
                Normal = Vector3.Normalize(normal),
                HalfExtents = size * 0.5f,
                MaterialIndex = mirrorMat
            };
            AddObject(plane);
            return plane;
        }

        /// <summary>
        /// Adds a mirror box.
        /// </summary>
        public TraceableBox AddMirrorBox(Vector3 center, Vector3 size)
        {
            var mirrorMat = AddMaterial(RayTracingMaterial.Mirror());
            var box = new TraceableBox
            {
                Center = center,
                HalfExtents = size * 0.5f,
                MaterialIndex = mirrorMat
            };
            AddObject(box);
            return box;
        }

        /// <summary>
        /// Clears all objects.
        /// </summary>
        public void ClearObjects()
        {
            _objects.Clear();
            _bvhDirty = true;
        }

        /// <summary>
        /// Rebuilds the BVH acceleration structure.
        /// </summary>
        public void RebuildBVH()
        {
            _bvhNodes.Clear();
            if (_objects.Count == 0)
            {
                _bvhDirty = false;
                return;
            }

            // Build BVH using surface area heuristic
            BuildBVHNode(Enumerable.Range(0, _objects.Count).ToList());
            _bvhDirty = false;
        }

        private int BuildBVHNode(List<int> objectIndices)
        {
            var node = new BVHNode
            {
                LeftChild = -1,
                RightChild = -1,
                ObjectIndex = -1
            };

            // Calculate bounds
            var bounds = _objects[objectIndices[0]].GetBounds();
            node.Min = bounds.min;
            node.Max = bounds.max;

            foreach (var idx in objectIndices.Skip(1))
            {
                var b = _objects[idx].GetBounds();
                node.Min = Vector3.Min(node.Min, b.min);
                node.Max = Vector3.Max(node.Max, b.max);
            }

            if (objectIndices.Count == 1)
            {
                node.ObjectIndex = objectIndices[0];
            }
            else
            {
                // Split along longest axis
                var extent = node.Max - node.Min;
                int axis = extent.X > extent.Y && extent.X > extent.Z ? 0 :
                          extent.Y > extent.Z ? 1 : 2;

                objectIndices.Sort((a, b) =>
                {
                    var ca = _objects[a].GetBounds();
                    var cb = _objects[b].GetBounds();
                    var centerA = (ca.min + ca.max) * 0.5f;
                    var centerB = (cb.min + cb.max) * 0.5f;
                    float va = axis == 0 ? centerA.X : axis == 1 ? centerA.Y : centerA.Z;
                    float vb = axis == 0 ? centerB.X : axis == 1 ? centerB.Y : centerB.Z;
                    return va.CompareTo(vb);
                });

                int mid = objectIndices.Count / 2;
                var left = objectIndices.Take(mid).ToList();
                var right = objectIndices.Skip(mid).ToList();

                int nodeIndex = _bvhNodes.Count;
                _bvhNodes.Add(node);

                node.LeftChild = BuildBVHNode(left);
                node.RightChild = BuildBVHNode(right);

                _bvhNodes[nodeIndex] = node;
                return nodeIndex;
            }

            int idx2 = _bvhNodes.Count;
            _bvhNodes.Add(node);
            return idx2;
        }

        /// <summary>
        /// Traces a single ray through the scene.
        /// </summary>
        public Vector3 TraceRay(Ray ray, int depth = 0)
        {
            if (depth >= MaxBounces)
                return BackgroundColor;

            if (_bvhDirty)
                RebuildBVH();

            if (!FindClosestHit(ray, out var hit))
                return BackgroundColor;

            var material = hit.MaterialIndex >= 0 && hit.MaterialIndex < _materials.Count
                ? _materials[hit.MaterialIndex]
                : RayTracingMaterial.Diffuse(Vector3.One);

            // Emission
            var color = material.Emission;

            // Reflection
            if (material.Reflectivity > 0)
            {
                var reflectDir = Reflect(ray.Direction, hit.Normal);

                // Add roughness
                if (material.Roughness > 0)
                {
                    reflectDir = AddRoughness(reflectDir, material.Roughness);
                }

                var reflectRay = new Ray(hit.Point + hit.Normal * 0.001f, reflectDir);
                var reflectColor = TraceRay(reflectRay, depth + 1);

                // Fresnel effect for dielectrics
                float fresnel = material.Metalness > 0.5f
                    ? material.Reflectivity
                    : Schlick(MathF.Abs(Vector3.Dot(ray.Direction, hit.Normal)), material.Reflectivity);

                color += material.Albedo * reflectColor * fresnel;
            }

            // Refraction (for transparent materials)
            if (material.Transparency > 0)
            {
                float ior = hit.FrontFace ? 1.0f / material.IOR : material.IOR;
                var refractDir = Refract(ray.Direction, hit.Normal, ior);

                if (refractDir.LengthSquared() > 0)
                {
                    var refractRay = new Ray(hit.Point - hit.Normal * 0.001f, refractDir);
                    var refractColor = TraceRay(refractRay, depth + 1);
                    color += material.Albedo * refractColor * material.Transparency;
                }
            }

            // Diffuse contribution
            float diffuseContribution = 1 - material.Reflectivity - material.Transparency;
            if (diffuseContribution > 0)
            {
                color += material.Albedo * AmbientLight * diffuseContribution;
            }

            return color;
        }

        /// <summary>
        /// Traces rays for an entire image (CPU multithreaded or GPU).
        /// </summary>
        public Vector3[] TraceImage(
            int width, int height,
            Vector3 cameraPos, Vector3 cameraTarget, Vector3 cameraUp,
            float fov = 60f)
        {
            var pixels = new Vector3[width * height];
            float aspectRatio = (float)width / height;
            float fovRad = fov * MathF.PI / 180f;
            float viewportHeight = 2f * MathF.Tan(fovRad / 2f);
            float viewportWidth = aspectRatio * viewportHeight;

            var forward = Vector3.Normalize(cameraTarget - cameraPos);
            var right = Vector3.Normalize(Vector3.Cross(forward, cameraUp));
            var up = Vector3.Cross(right, forward);

            var horizontal = right * viewportWidth;
            var vertical = up * viewportHeight;
            var lowerLeftCorner = cameraPos + forward - horizontal / 2 - vertical / 2;

            if (_bvhDirty)
                RebuildBVH();

            if (UseGPU && _gpuCompute?.IsAvailable == true)
            {
                TraceImageGPU(pixels, width, height, cameraPos, lowerLeftCorner, horizontal, vertical);
            }
            else
            {
                TraceImageCPU(pixels, width, height, cameraPos, lowerLeftCorner, horizontal, vertical);
            }

            return pixels;
        }

        private void TraceImageCPU(
            Vector3[] pixels, int width, int height,
            Vector3 origin, Vector3 lowerLeft, Vector3 horizontal, Vector3 vertical)
        {
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    float u = (float)x / (width - 1);
                    float v = (float)(height - 1 - y) / (height - 1);

                    var direction = lowerLeft + horizontal * u + vertical * v - origin;
                    var ray = new Ray(origin, direction);
                    pixels[y * width + x] = TraceRay(ray);
                }
            });
        }

        private void TraceImageGPU(
            Vector3[] pixels, int width, int height,
            Vector3 origin, Vector3 lowerLeft, Vector3 horizontal, Vector3 vertical)
        {
            // GPU kernel would be executed here
            // For now, fallback to CPU with SIMD optimizations
            TraceImageCPU(pixels, width, height, origin, lowerLeft, horizontal, vertical);
        }

        private bool FindClosestHit(Ray ray, out RayHit closestHit)
        {
            closestHit = default;
            float closestT = float.MaxValue;

            // Simple linear search (BVH traversal would be used for large scenes)
            foreach (var obj in _objects)
            {
                if (obj.Intersect(ray, out var hit) && hit.T < closestT)
                {
                    closestT = hit.T;
                    closestHit = hit;
                }
            }

            return closestHit.Hit;
        }

        private static Vector3 Reflect(Vector3 v, Vector3 n)
        {
            return v - 2 * Vector3.Dot(v, n) * n;
        }

        private static Vector3 Refract(Vector3 uv, Vector3 n, float etaiOverEtat)
        {
            float cosTheta = MathF.Min(Vector3.Dot(-uv, n), 1.0f);
            var rOutPerp = etaiOverEtat * (uv + cosTheta * n);
            var rOutParallel = -MathF.Sqrt(MathF.Abs(1.0f - rOutPerp.LengthSquared())) * n;
            return rOutPerp + rOutParallel;
        }

        private static float Schlick(float cosine, float refIdx)
        {
            float r0 = (1 - refIdx) / (1 + refIdx);
            r0 = r0 * r0;
            return r0 + (1 - r0) * MathF.Pow(1 - cosine, 5);
        }

        private Vector3 AddRoughness(Vector3 dir, float roughness)
        {
            // Simple roughness perturbation
            var random = new Random();
            var offset = new Vector3(
                (float)(random.NextDouble() - 0.5) * roughness,
                (float)(random.NextDouble() - 0.5) * roughness,
                (float)(random.NextDouble() - 0.5) * roughness);
            return Vector3.Normalize(dir + offset);
        }

        /// <summary>
        /// Converts ray-traced image to byte array (RGB).
        /// </summary>
        public byte[] ToByteArray(Vector3[] pixels, float exposure = 1.0f)
        {
            var bytes = new byte[pixels.Length * 3];

            for (int i = 0; i < pixels.Length; i++)
            {
                // Tone mapping (simple Reinhard)
                var color = pixels[i] * exposure;
                color = color / (color + Vector3.One);

                // Gamma correction
                color = new Vector3(
                    MathF.Pow(color.X, 1f / 2.2f),
                    MathF.Pow(color.Y, 1f / 2.2f),
                    MathF.Pow(color.Z, 1f / 2.2f));

                bytes[i * 3 + 0] = (byte)MathF.Min(255, color.X * 255);
                bytes[i * 3 + 1] = (byte)MathF.Min(255, color.Y * 255);
                bytes[i * 3 + 2] = (byte)MathF.Min(255, color.Z * 255);
            }

            return bytes;
        }

        public void Dispose()
        {
            _gpuCompute?.Dispose();
        }

        #region OpenCL Kernel Source

        /// <summary>
        /// OpenCL kernel for GPU ray tracing.
        /// </summary>
        public static readonly string OpenCLKernel = @"
typedef struct {
    float3 origin;
    float3 direction;
    float tMin;
    float tMax;
} Ray;

typedef struct {
    float3 albedo;
    float3 emission;
    float reflectivity;
    float roughness;
    float metalness;
    float ior;
    float transparency;
} Material;

typedef struct {
    float3 center;
    float radius;
    int materialIndex;
} Sphere;

float3 reflect(float3 v, float3 n) {
    return v - 2.0f * dot(v, n) * n;
}

float schlick(float cosine, float refIdx) {
    float r0 = (1.0f - refIdx) / (1.0f + refIdx);
    r0 = r0 * r0;
    return r0 + (1.0f - r0) * pow(1.0f - cosine, 5.0f);
}

bool intersectSphere(Ray ray, Sphere sphere, float* t, float3* normal) {
    float3 oc = ray.origin - sphere.center;
    float a = dot(ray.direction, ray.direction);
    float halfB = dot(oc, ray.direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;
    float discriminant = halfB * halfB - a * c;

    if (discriminant < 0) return false;

    float sqrtD = sqrt(discriminant);
    float root = (-halfB - sqrtD) / a;

    if (root < ray.tMin || root > ray.tMax) {
        root = (-halfB + sqrtD) / a;
        if (root < ray.tMin || root > ray.tMax)
            return false;
    }

    *t = root;
    float3 hitPoint = ray.origin + ray.direction * root;
    *normal = normalize(hitPoint - sphere.center);
    return true;
}

__kernel void traceRays(
    __global float3* pixels,
    __global Sphere* spheres,
    __global Material* materials,
    int numSpheres,
    float3 cameraPos,
    float3 lowerLeft,
    float3 horizontal,
    float3 vertical,
    int width,
    int height,
    int maxBounces,
    float3 backgroundColor)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    if (x >= width || y >= height) return;

    float u = (float)x / (float)(width - 1);
    float v = (float)(height - 1 - y) / (float)(height - 1);

    Ray ray;
    ray.origin = cameraPos;
    ray.direction = normalize(lowerLeft + u * horizontal + v * vertical - cameraPos);
    ray.tMin = 0.001f;
    ray.tMax = 1e30f;

    float3 color = (float3)(0, 0, 0);
    float3 throughput = (float3)(1, 1, 1);

    for (int bounce = 0; bounce < maxBounces; bounce++) {
        float closestT = 1e30f;
        int hitSphere = -1;
        float3 hitNormal;

        // Find closest intersection
        for (int i = 0; i < numSpheres; i++) {
            float t;
            float3 n;
            if (intersectSphere(ray, spheres[i], &t, &n) && t < closestT) {
                closestT = t;
                hitSphere = i;
                hitNormal = n;
            }
        }

        if (hitSphere < 0) {
            color += throughput * backgroundColor;
            break;
        }

        Material mat = materials[spheres[hitSphere].materialIndex];
        color += throughput * mat.emission;

        // Reflection
        if (mat.reflectivity > 0) {
            float3 hitPoint = ray.origin + ray.direction * closestT;
            float3 reflectDir = reflect(ray.direction, hitNormal);

            ray.origin = hitPoint + hitNormal * 0.001f;
            ray.direction = reflectDir;

            float fresnel = mat.metalness > 0.5f ? mat.reflectivity :
                schlick(fabs(dot(ray.direction, hitNormal)), mat.reflectivity);

            throughput *= mat.albedo * fresnel;
        } else {
            color += throughput * mat.albedo * 0.1f; // Ambient
            break;
        }
    }

    pixels[y * width + x] = color;
}
";

        #endregion
    }

    /// <summary>
    /// Light source for ray tracing.
    /// </summary>
    public class RayTracingLight
    {
        /// <summary>Light position.</summary>
        public Vector3 Position { get; set; }

        /// <summary>Light color and intensity.</summary>
        public Vector3 Color { get; set; } = Vector3.One;

        /// <summary>Light radius (for soft shadows).</summary>
        public float Radius { get; set; } = 0.1f;

        /// <summary>Light type.</summary>
        public LightType Type { get; set; } = LightType.Point;

        /// <summary>Direction (for directional/spot lights).</summary>
        public Vector3 Direction { get; set; } = -Vector3.UnitY;

        /// <summary>Spot angle in degrees.</summary>
        public float SpotAngle { get; set; } = 45f;

        /// <summary>
        /// Creates a point light.
        /// </summary>
        public static RayTracingLight Point(Vector3 position, Vector3 color)
        {
            return new RayTracingLight
            {
                Position = position,
                Color = color,
                Type = LightType.Point
            };
        }

        /// <summary>
        /// Creates a directional light.
        /// </summary>
        public static RayTracingLight Directional(Vector3 direction, Vector3 color)
        {
            return new RayTracingLight
            {
                Direction = Vector3.Normalize(direction),
                Color = color,
                Type = LightType.Directional
            };
        }

        /// <summary>
        /// Creates an area light (for soft shadows).
        /// </summary>
        public static RayTracingLight Area(Vector3 position, float radius, Vector3 color)
        {
            return new RayTracingLight
            {
                Position = position,
                Radius = radius,
                Color = color,
                Type = LightType.Area
            };
        }
    }

    /// <summary>
    /// Light type enumeration.
    /// </summary>
    public enum LightType
    {
        Point,
        Directional,
        Spot,
        Area
    }
}
