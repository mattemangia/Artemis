namespace Artemis.Compute
{
    /// <summary>
    /// OpenCL C source of the particle solver. It is a line-by-line port of the CPU
    /// solver in <c>MassiveParticleSystem.Cpu.cs</c>; keep the two in step.
    /// </summary>
    internal static class ParticleKernels
    {
        /// <summary>Work-group size of the prefix-sum kernels (fixed by the source).</summary>
        public const int ScanGroupSize = 256;

        /// <summary>Elements handled by one prefix-sum work group.</summary>
        public const int ScanBlock = ScanGroupSize * 4;

        public const string Source = """
#define FLAG_ALIVE 1
#define FLAG_COLLIDES 2
#define FLAG_GRAVITY 4
#define FLAG_FROZEN 8
#define FLAG_SLEEPING 16
#define WAKE_UP 1
#define WAKE_RELEASE 2
#define SOFT_TAG 0x40000000

inline int cell_coord(float v, float inv) { return (int)floor(v * inv); }

/* Rows of cells (fixed y and z) are hashed; the cells of a row occupy consecutive
   buckets, so the three cells x-1, x, x+1 are one contiguous run. */
inline uint hash_row(int cy, int cz) { return ((uint)cy * 19349663u) ^ ((uint)cz * 83492791u); }
inline int cell_key(int cx, int cy, int cz) { return (cx & 1023) | ((cy & 1023) << 10) | ((cz & 1023) << 20); }

/* ---------------------------------------------------------------- predict */

__kernel void predict(int n, float dt, float gxdt, float gydt, float gzdt, float damping,
    __global float* px, __global float* py, __global float* pz,
    __global float* qx, __global float* qy, __global float* qz,
    __global float* vx, __global float* vy, __global float* vz,
    __global const float* invMass, __global uchar* flags, __global uchar* rest, __global int* wake,
    __global float* dpx, __global float* dpy, __global float* dpz,
    __global float* stx, __global float* sty, __global float* stz,
    __global int* counters)
{
    __local int localAwake;
    if (get_local_id(0) == 0) localAwake = 0;
    barrier(CLK_LOCAL_MEM_FENCE);

    int i = get_global_id(0);
    int awake = 0;
    if (i < n)
    {
        uchar f = flags[i];
        int w = wake[i];
        if (w != 0)
        {
            wake[i] = 0;
            f &= (uchar)~FLAG_SLEEPING;
            if (w == WAKE_RELEASE) f &= (uchar)~FLAG_FROZEN;
            flags[i] = f;
            rest[i] = 0;
        }

        float x = px[i], y = py[i], z = pz[i];
        qx[i] = x; qy[i] = y; qz[i] = z;
        dpx[i] = 0; dpy[i] = 0; dpz[i] = 0;
        stx[i] = 0; sty[i] = 0; stz[i] = 0;

        if ((f & FLAG_ALIVE) != 0 && (f & (FLAG_FROZEN | FLAG_SLEEPING)) == 0 && invMass[i] > 0)
        {
            if ((f & FLAG_COLLIDES) != 0) awake = 1;

            float ux = vx[i], uy = vy[i], uz = vz[i];
            if ((f & FLAG_GRAVITY) != 0) { ux += gxdt; uy += gydt; uz += gzdt; }
            ux *= damping; uy *= damping; uz *= damping;
            vx[i] = ux; vy[i] = uy; vz[i] = uz;
            px[i] = x + ux * dt; py[i] = y + uy * dt; pz[i] = z + uz * dt;
        }
    }

    if (awake) atomic_inc(&localAwake);
    barrier(CLK_LOCAL_MEM_FENCE);
    if (get_local_id(0) == 0 && localAwake != 0) atomic_add(&counters[0], localAwake);
}

/* ---------------------------------------------------------------- grid */

__kernel void assign_cells(int n, float invCell, int mask,
    __global const float* px, __global const float* py, __global const float* pz,
    __global const uchar* flags, __global int* cell, __global int* bucketCount)
{
    int i = get_global_id(0);
    if (i >= n) return;

    if ((flags[i] & (FLAG_ALIVE | FLAG_COLLIDES)) != (FLAG_ALIVE | FLAG_COLLIDES))
    {
        cell[i] = -1;
        return;
    }

    int cx = cell_coord(px[i], invCell);
    int cy = cell_coord(py[i], invCell);
    int cz = cell_coord(pz[i], invCell);
    int b = (int)((hash_row(cy, cz) + (uint)cx) & (uint)mask);
    cell[i] = b;
    atomic_inc(&bucketCount[b]);
}

/* Exclusive prefix sum of 1024 elements per work group (256 items x 4). */
__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void scan_block(__global const int* input, __global int* output, __global int* blockSums, int n)
{
    __local int temp[256];
    int lid = get_local_id(0);
    int base = get_group_id(0) * 1024 + lid * 4;

    int v0 = base < n ? input[base] : 0;
    int v1 = base + 1 < n ? input[base + 1] : 0;
    int v2 = base + 2 < n ? input[base + 2] : 0;
    int v3 = base + 3 < n ? input[base + 3] : 0;
    int s1 = v0 + v1, s2 = s1 + v2, s3 = s2 + v3;

    temp[lid] = s3;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int offset = 1; offset < 256; offset <<= 1)
    {
        int t = lid >= offset ? temp[lid - offset] : 0;
        barrier(CLK_LOCAL_MEM_FENCE);
        temp[lid] += t;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    int exclusive = temp[lid] - s3;
    if (base < n) output[base] = exclusive;
    if (base + 1 < n) output[base + 1] = exclusive + v0;
    if (base + 2 < n) output[base + 2] = exclusive + s1;
    if (base + 3 < n) output[base + 3] = exclusive + s2;
    if (lid == 255) blockSums[get_group_id(0)] = temp[255];
}

__kernel void add_offsets(__global int* data, __global const int* scannedSums, int n)
{
    int i = get_global_id(0);
    if (i < n) data[i] += scannedSums[i / 1024];
}

/* start[tableSize .. tableSize + 3] = total, so that a run of three buckets can be
   read as start[b] .. start[b + 3] even at the end of the table. */
__kernel void finalize_starts(__global int* start, __global const int* count, int tableSize)
{
    if (get_global_id(0) != 0) return;
    int total = start[tableSize - 1] + count[tableSize - 1];
    start[tableSize] = total;
    start[tableSize + 1] = total;
    start[tableSize + 2] = total;
    start[tableSize + 3] = total;
}

__kernel void scatter(int n, __global const int* cell, __global int* cursor, __global int* sortedIndex)
{
    int i = get_global_id(0);
    if (i >= n) return;
    int b = cell[i];
    if (b >= 0) sortedIndex[atomic_inc(&cursor[b])] = i;
}

/* Atomics scatter in an arbitrary order; sorting each bucket by particle index
   makes the summation order, and with it the result, reproducible. */
__kernel void sort_buckets(int tableSize, __global const int* start, __global int* sortedIndex)
{
    int b = get_global_id(0);
    if (b >= tableSize) return;
    int from = start[b], to = start[b + 1];
    for (int a = from + 1; a < to; a++)
    {
        int value = sortedIndex[a];
        int c = a - 1;
        while (c >= from && sortedIndex[c] > value)
        {
            sortedIndex[c + 1] = sortedIndex[c];
            c--;
        }
        sortedIndex[c + 1] = value;
    }
}

__kernel void gather(int mask, float invCell, __global const int* start, __global const int* sortedIndex,
    __global const float* px, __global const float* py, __global const float* pz,
    __global const float* qx, __global const float* qy, __global const float* qz,
    __global const float* radius, __global const float* invMass, __global const uchar* flags, __global const int* body,
    __global float* sx, __global float* sy, __global float* sz, __global float* sr, __global float* sw,
    __global float* smx, __global float* smy, __global float* smz, __global int* skey, __global int* sbody)
{
    int s = get_global_id(0);
    if (s >= start[mask + 1]) return;

    int i = sortedIndex[s];
    uchar f = flags[i];
    float x = px[i], y = py[i], z = pz[i];
    sx[s] = x; sy[s] = y; sz[s] = z;
    sr[s] = radius[i];
    skey[s] = cell_key(cell_coord(x, invCell), cell_coord(y, invCell), cell_coord(z, invCell));
    sw[s] = (f & (FLAG_FROZEN | FLAG_SLEEPING)) == 0 ? invMass[i] : ((f & FLAG_FROZEN) != 0 ? -1.0f : -2.0f);
    sbody[s] = body[i];
    smx[s] = x - qx[i]; smy[s] = y - qy[i]; smz[s] = z - qz[i];
}

/* ---------------------------------------------------------------- contacts */

__kernel void solve(int mask, float invCell, float friction, float relaxation, float movingThresholdSq,
    float breakDistance, int canBreak, float stackingBias, float upX, float upY, float upZ,
    __global const int* start, __global const int* sortedIndex,
    __global const float* sx, __global const float* sy, __global const float* sz,
    __global const float* sr, __global const float* sw,
    __global const float* smx, __global const float* smy, __global const float* smz,
    __global const int* skey, __global const int* sbody,
    __global int* wake,
    __global float* dpx, __global float* dpy, __global float* dpz,
    __global float* stx, __global float* sty, __global float* stz)
{
    int s = get_global_id(0);
    if (s >= start[mask + 1]) return;

    float wi = sw[s];
    if (wi <= 0) return;

    int selfBody = sbody[s];
    int selfIndex = sortedIndex[s];
    float px = sx[s], py = sy[s], pz = sz[s];
    float ri = sr[s];
    float mx = smx[s], my = smy[s], mz = smz[s];
    int moving = mx * mx + my * my + mz * mz > movingThresholdSq;

    int cx = cell_coord(px, invCell);
    int cy = cell_coord(py, invCell);
    int cz = cell_coord(pz, invCell);

    float dx_ = 0, dy_ = 0, dz_ = 0;
    float tx_ = 0, ty_ = 0, tz_ = 0;
    int contacts = 0;
    float nxx = 0, nxy = 0, nxz = 0, nyy = 0, nyz = 0, nzz = 0;
    float fx = 0, fy = 0, fz = 0;

    for (int oz = -1; oz <= 1; oz++)
    {
        for (int oy = -1; oy <= 1; oy++)
        {
            int rowKey = (((cy + oy) & 1023) << 10) | (((cz + oz) & 1023) << 20);
            int first = (int)((hash_row(cy + oy, cz + oz) + (uint)(cx - 1)) & (uint)mask);
            int pieces = first <= mask - 2 ? 1 : 3;

            for (int piece = 0; piece < pieces; piece++)
            {
                int from, to;
                if (pieces == 1) { from = start[first]; to = start[first + 3]; }
                else { int b = (first + piece) & mask; from = start[b]; to = start[b + 1]; }

                for (int t = from; t < to; t++)
                {
                    if (selfBody >= 0 && sbody[t] == selfBody) continue;
                    int key = skey[t];
                    if (t == s || ((key ^ rowKey) & ~1023) != 0 || ((key - cx + 1) & 1023) > 2) continue;

                    float dx = px - sx[t], dy = py - sy[t], dz = pz - sz[t];
                    float minDist = ri + sr[t];
                    float distSq = dx * dx + dy * dy + dz * dz;
                    float wj = sw[t];

                    if (wj == -2.0f && moving && distSq < minDist * minDist * 1.21f)
                        atomic_max(&wake[sortedIndex[t]], WAKE_UP);

                    if (distSq >= minDist * minDist) continue;

                    float nx, ny, nz, dist;
                    if (distSq > 1e-12f)
                    {
                        dist = sqrt(distSq);
                        float inv = 1.0f / dist;
                        nx = dx * inv; ny = dy * inv; nz = dz * inv;
                    }
                    else
                    {
                        dist = 0;
                        nx = 0; ny = selfIndex > sortedIndex[t] ? 1.0f : -1.0f; nz = 0;
                    }

                    float rx = mx, ry = my, rz = mz;
                    if (wj > 0) { rx -= smx[t]; ry -= smy[t]; rz -= smz[t]; }
                    float rn = rx * nx + ry * ny + rz * nz;

                    if (wj == -1.0f && canBreak && -rn > breakDistance)
                        atomic_max(&wake[sortedIndex[t]], WAKE_RELEASE);

                    float penetration = minDist - dist;
                    float wOther = wj > 0 ? wj : 0;
                    float share = wi / (wi + wOther);
                    if (wOther > 0 && stackingBias > 0)
                    {
                        float up = nx * upX + ny * upY + nz * upZ;
                        share = clamp(share * (1.0f + stackingBias * up), 0.0f, 1.0f);
                    }

                    float magnitude = penetration * share;
                    dx_ += nx * magnitude; dy_ += ny * magnitude; dz_ += nz * magnitude;
                    nxx += nx * nx; nxy += nx * ny; nxz += nx * nz;
                    nyy += ny * ny; nyz += ny * nz; nzz += nz * nz;
                    contacts++;

                    float approach = rn < 0 ? -rn : 0;
                    float preexisting = penetration - approach;
                    if (preexisting > 0)
                    {
                        tx_ += nx * preexisting * share;
                        ty_ += ny * preexisting * share;
                        tz_ += nz * preexisting * share;
                    }

                    if (friction > 0)
                    {
                        float tx = rx - rn * nx, ty = ry - rn * ny, tz = rz - rn * nz;
                        float slideSq = tx * tx + ty * ty + tz * tz;
                        if (slideSq > 1e-16f)
                        {
                            float slide = sqrt(slideSq);
                            float limit = friction * penetration;
                            float scale = (slide <= limit ? 1.0f : limit / slide) * share;
                            fx -= tx * scale; fy -= ty * scale; fz -= tz * scale;
                        }
                    }
                }
            }
        }
    }

    if (contacts > 0)
    {
        float lengthSq = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;
        float projected =
            dx_ * (nxx * dx_ + nxy * dy_ + nxz * dz_) +
            dy_ * (nxy * dx_ + nyy * dy_ + nyz * dz_) +
            dz_ * (nxz * dx_ + nyz * dy_ + nzz * dz_);
        float alpha = projected > 1e-20f ? fmin(lengthSq / projected, 1.5f) : 0.0f;
        float scale = relaxation * alpha;

        float frictionScale = fmin(alpha, 1.0f);
        fx *= frictionScale; fy *= frictionScale; fz *= frictionScale;
        float frictionSq = fx * fx + fy * fy + fz * fz;
        float movedSq = mx * mx + my * my + mz * mz;
        if (frictionSq > movedSq && frictionSq > 0)
        {
            float c = sqrt(movedSq / frictionSq);
            fx *= c; fy *= c; fz *= c;
        }

        dpx[selfIndex] = dx_ * scale + fx;
        dpy[selfIndex] = dy_ * scale + fy;
        dpz[selfIndex] = dz_ * scale + fz;
        stx[selfIndex] = tx_ * scale;
        sty[selfIndex] = ty_ * scale;
        stz[selfIndex] = tz_ * scale;
    }
}

__kernel void apply_contacts(int mask, __global const int* start, __global const int* sortedIndex, __global const float* sw,
    __global float* px, __global float* py, __global float* pz,
    __global const float* dpx, __global const float* dpy, __global const float* dpz)
{
    int s = get_global_id(0);
    if (s >= start[mask + 1] || sw[s] <= 0) return;
    int i = sortedIndex[s];
    px[i] += dpx[i]; py[i] += dpy[i]; pz[i] += dpz[i];
}

/* ---------------------------------------------------------------- soft bodies */

inline float4 quat_mul(float4 a, float4 b)
{
    return (float4)(
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
}

inline float3 quat_rotate(float4 q, float3 v)
{
    float3 u = q.xyz;
    float3 t = 2.0f * cross(u, v);
    return v + q.w * t + cross(u, t);
}

/* One cluster per work item: best-fit rotation of the cluster (iterative method of
   Mueller et al. 2016, warm-started from the previous rotation) and the goal
   position of every member. */
__kernel void soft_fit(int clusterCount,
    __global const int* clusterStart, __global const int* members,
    __global const float* restX, __global const float* restY, __global const float* restZ,
    __global float4* rotation,
    __global const float* px, __global const float* py, __global const float* pz, __global const float* invMass,
    __global const uchar* flags,
    __global float* goalX, __global float* goalY, __global float* goalZ)
{
    int c = get_global_id(0);
    if (c >= clusterCount) return;

    int from = clusterStart[c], to = clusterStart[c + 1];
    float total = 0;
    float3 center = (float3)(0);
    for (int k = from; k < to; k++)
    {
        int i = members[k];
        float m = (flags[i] & FLAG_ALIVE) != 0 && invMass[i] > 0 ? 1.0f / invMass[i] : 0.0f;
        center += m * (float3)(px[i], py[i], pz[i]);
        total += m;
    }
    if (total <= 0) return;
    center /= total;

    float3 a0 = (float3)(0), a1 = (float3)(0), a2 = (float3)(0);
    for (int k = from; k < to; k++)
    {
        int i = members[k];
        float m = (flags[i] & FLAG_ALIVE) != 0 && invMass[i] > 0 ? 1.0f / invMass[i] : 0.0f;
        float3 d = m * ((float3)(px[i], py[i], pz[i]) - center);
        a0 += d * restX[k];
        a1 += d * restY[k];
        a2 += d * restZ[k];
    }

    float4 q = rotation[c];
    for (int iteration = 0; iteration < 20; iteration++)
    {
        float3 r0 = quat_rotate(q, (float3)(1, 0, 0));
        float3 r1 = quat_rotate(q, (float3)(0, 1, 0));
        float3 r2 = quat_rotate(q, (float3)(0, 0, 1));
        float3 omega = cross(r0, a0) + cross(r1, a1) + cross(r2, a2);
        omega /= fabs(dot(r0, a0) + dot(r1, a1) + dot(r2, a2)) + 1e-12f;
        float angle = length(omega);
        if (angle < 1e-7f) break;
        float3 axis = omega / angle;
        float s = sin(angle * 0.5f);
        q = normalize(quat_mul((float4)(axis * s, cos(angle * 0.5f)), q));
    }
    rotation[c] = q;

    for (int k = from; k < to; k++)
    {
        float3 g = center + quat_rotate(q, (float3)(restX[k], restY[k], restZ[k]));
        goalX[k] = g.x; goalY[k] = g.y; goalZ[k] = g.z;
    }
}

/* One particle per work item: move towards the mean of its clusters' goals. */
__kernel void soft_apply(int n, __global const int* slotStart, __global const int* slots, __global const int* body,
    __global const float* alpha,
    __global const float* goalX, __global const float* goalY, __global const float* goalZ,
    __global float* px, __global float* py, __global float* pz)
{
    int i = get_global_id(0);
    if (i >= n) return;
    int from = slotStart[i], to = slotStart[i + 1];
    if (to <= from) return;

    float3 sum = (float3)(0);
    for (int k = from; k < to; k++)
    {
        int slot = slots[k];
        sum += (float3)(goalX[slot], goalY[slot], goalZ[slot]);
    }
    float3 p = (float3)(px[i], py[i], pz[i]);
    float3 d = sum / (float)(to - from) - p;
    float3 c = alpha[body[i] & ~SOFT_TAG] * d;
    px[i] = p.x + c.x; py[i] = p.y + c.y; pz[i] = p.z + c.z;
}

/* Internal damping (Mueller et al. 2007): velocities are pulled towards the rigid
   motion of their body. Step 1: sums over chunks of 256 members, relative to the
   body's first particle. */
__kernel void soft_body_sum(int chunkCount, __global const int* chunkBody, __global const int* chunkFrom, __global const int* chunkTo,
    __global const int* bodyStart, __global const int* bodyMembers,
    __global const float* px, __global const float* py, __global const float* pz,
    __global const float* vx, __global const float* vy, __global const float* vz,
    __global const float* invMass, __global const uchar* flags, __global float* partial)
{
    int c = get_global_id(0);
    if (c >= chunkCount) return;
    int ref = bodyMembers[bodyStart[chunkBody[c]]];
    float3 o = (float3)(px[ref], py[ref], pz[ref]);

    float m = 0;
    float3 mr = (float3)(0), p = (float3)(0), l = (float3)(0);
    float sxx = 0, sxy = 0, sxz = 0, syy = 0, syz = 0, szz = 0;
    for (int k = chunkFrom[c]; k < chunkTo[c]; k++)
    {
        int i = bodyMembers[k];
        float w = (flags[i] & FLAG_ALIVE) != 0 && invMass[i] > 0 ? 1.0f / invMass[i] : 0.0f;
        float3 r = (float3)(px[i], py[i], pz[i]) - o;
        float3 v = (float3)(vx[i], vy[i], vz[i]);
        m += w; mr += w * r; p += w * v; l += w * cross(r, v);
        sxx += w * r.x * r.x; sxy += w * r.x * r.y; sxz += w * r.x * r.z;
        syy += w * r.y * r.y; syz += w * r.y * r.z; szz += w * r.z * r.z;
    }

    __global float* out = partial + 16 * c;
    out[0] = m; out[1] = mr.x; out[2] = mr.y; out[3] = mr.z;
    out[4] = p.x; out[5] = p.y; out[6] = p.z;
    out[7] = l.x; out[8] = l.y; out[9] = l.z;
    out[10] = sxx; out[11] = sxy; out[12] = sxz; out[13] = syy; out[14] = syz; out[15] = szz;
}

/* Step 2: per body, centre-of-mass velocity, centre and angular velocity. */
__kernel void soft_body_rigid(int bodyCount, __global const int* bodyChunkStart, __global const float* partial,
    __global const int* bodyStart, __global const int* bodyMembers,
    __global const float* px, __global const float* py, __global const float* pz,
    __global const float* beta, __global float* rigid)
{
    int b = get_global_id(0);
    if (b >= bodyCount) return;

    __global float* out = rigid + 9 * b;
    for (int k = 0; k < 9; k++) out[k] = 0;
    if (beta[b] <= 0 || bodyChunkStart[b + 1] <= bodyChunkStart[b]) return;

    float s[16];
    for (int k = 0; k < 16; k++) s[k] = 0;
    for (int c = bodyChunkStart[b]; c < bodyChunkStart[b + 1]; c++)
        for (int k = 0; k < 16; k++) s[k] += partial[16 * c + k];

    float m = s[0];
    if (m <= 0) return;
    float3 d = (float3)(s[1], s[2], s[3]) / m;
    float3 vcm = (float3)(s[4], s[5], s[6]) / m;
    float3 l = (float3)(s[7], s[8], s[9]) - m * cross(d, vcm);
    float sxx = s[10] - m * d.x * d.x, sxy = s[11] - m * d.x * d.y, sxz = s[12] - m * d.x * d.z;
    float syy = s[13] - m * d.y * d.y, syz = s[14] - m * d.y * d.z, szz = s[15] - m * d.z * d.z;

    float trace = sxx + syy + szz;
    float eps = 1e-6f * trace + 1e-12f;
    float ixx = syy + szz + eps, iyy = sxx + szz + eps, izz = sxx + syy + eps;
    float ixy = -sxy, ixz = -sxz, iyz = -syz;
    float c00 = iyy * izz - iyz * iyz, c01 = ixz * iyz - ixy * izz, c02 = ixy * iyz - ixz * iyy;
    float det = ixx * c00 + ixy * c01 + ixz * c02;

    int ref = bodyMembers[bodyStart[b]];
    out[0] = vcm.x; out[1] = vcm.y; out[2] = vcm.z;
    out[3] = px[ref] + d.x; out[4] = py[ref] + d.y; out[5] = pz[ref] + d.z;
    if (fabs(det) < 1e-30f) return;

    float c11 = ixx * izz - ixz * ixz, c12 = ixy * ixz - ixx * iyz, c22 = ixx * iyy - ixy * ixy;
    float inv = 1.0f / det;
    out[6] = (c00 * l.x + c01 * l.y + c02 * l.z) * inv;
    out[7] = (c01 * l.x + c11 * l.y + c12 * l.z) * inv;
    out[8] = (c02 * l.x + c12 * l.y + c22 * l.z) * inv;
}

/* Step 3: pull each soft particle's velocity towards its body's rigid motion. */
__kernel void soft_body_damp(int n, __global const int* body, __global const float* beta, __global const float* rigid,
    __global const uchar* flags,
    __global const float* px, __global const float* py, __global const float* pz,
    __global float* vx, __global float* vy, __global float* vz)
{
    int i = get_global_id(0);
    if (i >= n) return;
    int tag = body[i];
    if (tag < 0 || (tag & SOFT_TAG) == 0 || (flags[i] & (FLAG_ALIVE | FLAG_FROZEN | FLAG_SLEEPING)) != FLAG_ALIVE) return;
    int b = tag & ~SOFT_TAG;
    float k = beta[b];
    if (k <= 0) return;

    __global const float* r = rigid + 9 * b;
    float3 rel = (float3)(px[i] - r[3], py[i] - r[4], pz[i] - r[5]);
    float3 target = (float3)(r[0], r[1], r[2]) + cross((float3)(r[6], r[7], r[8]), rel);
    vx[i] += k * (target.x - vx[i]);
    vy[i] += k * (target.y - vy[i]);
    vz[i] += k * (target.z - vz[i]);
}

/* ---------------------------------------------------------------- finish */

inline void slide(float* a, float* b, float moveA, float moveB, float depth, float friction)
{
    if (friction <= 0) return;
    float slideSq = moveA * moveA + moveB * moveB;
    if (slideSq < 1e-16f) return;
    float s = sqrt(slideSq);
    float limit = friction * depth;
    float scale = s <= limit ? 1.0f : limit / s;
    *a -= moveA * scale;
    *b -= moveB * scale;
}

__kernel void finish(int n, float dt, float minX, float minY, float minZ, float maxX, float maxY, float maxZ,
    float restitution, float friction, int allowSleeping, float sleepSpeedSq, int sleepSubSteps, int solved,
    __global float* px, __global float* py, __global float* pz,
    __global const float* qx, __global const float* qy, __global const float* qz,
    __global float* vx, __global float* vy, __global float* vz,
    __global const float* invMass, __global const float* radius, __global float* lifetime,
    __global uchar* flags, __global uchar* rest, __global const int* body,
    __global const float* dpx, __global const float* dpy, __global const float* dpz,
    __global const float* stx, __global const float* sty, __global const float* stz,
    __global int* counters)
{
    __local int localActive, localAlive;
    if (get_local_id(0) == 0) { localActive = 0; localAlive = 0; }
    barrier(CLK_LOCAL_MEM_FENCE);

    int i = get_global_id(0);
    int isActive = 0, isAlive = 0;
    const float invDt = 1.0f / dt;
    const float bounceThreshold = 1.0f;

    if (i < n && (flags[i] & FLAG_ALIVE) != 0)
    {
        uchar f = flags[i];
        isAlive = 1;

        if (lifetime[i] != FLT_MAX)
        {
            lifetime[i] -= dt;
            if (lifetime[i] <= 0)
            {
                flags[i] = (uchar)(f & ~FLAG_ALIVE);
                isAlive = 0;
            }
        }

        if (isAlive && (f & (FLAG_FROZEN | FLAG_SLEEPING)) == 0 && invMass[i] > 0)
        {
            isActive = 1;
            float r = radius[i];
            float x = px[i], y = py[i], z = pz[i];
            float moveX = x - qx[i], moveY = y - qy[i], moveZ = z - qz[i];
            int hit = 0;
            float depth;

            if ((depth = minX + r - x) > 0) { x = minX + r; hit |= 1; slide(&y, &z, moveY, moveZ, depth, friction); }
            else if ((depth = x - (maxX - r)) > 0) { x = maxX - r; hit |= 1; slide(&y, &z, moveY, moveZ, depth, friction); }
            if ((depth = minY + r - y) > 0) { y = minY + r; hit |= 2; slide(&x, &z, moveX, moveZ, depth, friction); }
            else if ((depth = y - (maxY - r)) > 0) { y = maxY - r; hit |= 2; slide(&x, &z, moveX, moveZ, depth, friction); }
            if ((depth = minZ + r - z) > 0) { z = minZ + r; hit |= 4; slide(&x, &y, moveX, moveY, depth, friction); }
            else if ((depth = z - (maxZ - r)) > 0) { z = maxZ - r; hit |= 4; slide(&x, &y, moveX, moveY, depth, friction); }

            px[i] = x; py[i] = y; pz[i] = z;

            int corrected = hit != 0 || body[i] >= 0 ||
                (solved && (f & FLAG_COLLIDES) != 0 && (dpx[i] != 0 || dpy[i] != 0 || dpz[i] != 0));

            if (!corrected)
            {
                if (allowSleeping) rest[i] = 0;
            }
            else
            {
                float ux = (x - qx[i]) * invDt, uy = (y - qy[i]) * invDt, uz = (z - qz[i]) * invDt;
                if (solved && (f & FLAG_COLLIDES) != 0)
                {
                    float ax = ux - stx[i] * invDt, ay = uy - sty[i] * invDt, az = uz - stz[i] * invDt;
                    if (ax * ax + ay * ay + az * az < ux * ux + uy * uy + uz * uz) { ux = ax; uy = ay; uz = az; }
                }

                if (hit != 0 && restitution > 0)
                {
                    if ((hit & 1) != 0 && fabs(vx[i]) > bounceThreshold) ux = -vx[i] * restitution;
                    if ((hit & 2) != 0 && fabs(vy[i]) > bounceThreshold) uy = -vy[i] * restitution;
                    if ((hit & 4) != 0 && fabs(vz[i]) > bounceThreshold) uz = -vz[i] * restitution;
                }

                vx[i] = ux; vy[i] = uy; vz[i] = uz;

                if (allowSleeping && body[i] < 0)
                {
                    if (ux * ux + uy * uy + uz * uz < sleepSpeedSq)
                    {
                        int r2 = rest[i] + 1;
                        if (r2 >= sleepSubSteps)
                        {
                            flags[i] = (uchar)(f | FLAG_SLEEPING);
                            vx[i] = 0; vy[i] = 0; vz[i] = 0;
                            r2 = 0;
                        }
                        rest[i] = (uchar)r2;
                    }
                    else
                    {
                        rest[i] = 0;
                    }
                }
            }
        }
    }

    if (isActive) atomic_inc(&localActive);
    if (isAlive) atomic_inc(&localAlive);
    barrier(CLK_LOCAL_MEM_FENCE);
    if (get_local_id(0) == 0)
    {
        if (localActive != 0) atomic_add(&counters[1], localActive);
        if (localAlive != 0) atomic_add(&counters[2], localAlive);
    }
}
""";
    }
}
