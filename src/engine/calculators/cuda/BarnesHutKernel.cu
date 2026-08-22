#include "BarnesHutKernel.cuh"

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

namespace galaxy_cuda
{

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            throw std::runtime_error(                                          \
                std::string("CUDA error: ") + cudaGetErrorString(err) +        \
                " (" + __FILE__ + ":" + std::to_string(__LINE__) + ")");       \
        }                                                                      \
    } while (0)

constexpr int EMPTY     = -1;
constexpr int STACK_CAP = 64;

struct BarnesHutDeviceBuffers
{
    int capacityParticles = 0;
    int capacityNodes     = 0;

    float3* dPos          = nullptr;
    float*  dPartMass     = nullptr;
    float3* dAcc          = nullptr;

    float3* dCom          = nullptr;
    float*  dNodeMass     = nullptr;
    float*  dHalfSize     = nullptr;
    int*    dChildren     = nullptr;
    int*    dParticle     = nullptr;
    int*    dHasChildren  = nullptr;
    int*    dDepth        = nullptr;

    void ensure(int nParticles, int nodeCount)
    {
        if (nParticles > capacityParticles)
        {
            freeParticles();
            capacityParticles = std::max(nParticles, capacityParticles * 2);
            if (capacityParticles < nParticles) capacityParticles = nParticles;
            CUDA_CHECK(cudaMalloc(&dPos,      capacityParticles * sizeof(float3)));
            CUDA_CHECK(cudaMalloc(&dPartMass, capacityParticles * sizeof(float)));
            CUDA_CHECK(cudaMalloc(&dAcc,      capacityParticles * sizeof(float3)));
        }
        if (nodeCount > capacityNodes)
        {
            freeNodes();
            capacityNodes = std::max(nodeCount, capacityNodes * 2);
            if (capacityNodes < nodeCount) capacityNodes = nodeCount;
            CUDA_CHECK(cudaMalloc(&dCom,         capacityNodes * sizeof(float3)));
            CUDA_CHECK(cudaMalloc(&dNodeMass,    capacityNodes * sizeof(float)));
            CUDA_CHECK(cudaMalloc(&dHalfSize,    capacityNodes * sizeof(float)));
            CUDA_CHECK(cudaMalloc(&dChildren,    capacityNodes * 8 * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&dParticle,    capacityNodes * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&dHasChildren, capacityNodes * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&dDepth,       capacityNodes * sizeof(int)));
        }
    }

    void freeParticles()
    {
        if (dPos)      { cudaFree(dPos);      dPos = nullptr; }
        if (dPartMass) { cudaFree(dPartMass); dPartMass = nullptr; }
        if (dAcc)      { cudaFree(dAcc);      dAcc = nullptr; }
        capacityParticles = 0;
    }

    void freeNodes()
    {
        if (dCom)         { cudaFree(dCom);         dCom = nullptr; }
        if (dNodeMass)    { cudaFree(dNodeMass);    dNodeMass = nullptr; }
        if (dHalfSize)    { cudaFree(dHalfSize);    dHalfSize = nullptr; }
        if (dChildren)    { cudaFree(dChildren);    dChildren = nullptr; }
        if (dParticle)    { cudaFree(dParticle);    dParticle = nullptr; }
        if (dHasChildren) { cudaFree(dHasChildren); dHasChildren = nullptr; }
        if (dDepth)       { cudaFree(dDepth);       dDepth = nullptr; }
        capacityNodes = 0;
    }

    void release()
    {
        freeParticles();
        freeNodes();
    }
};

static BarnesHutDeviceBuffers gBuffers;

void releaseBarnesHutCudaBuffers()
{
    gBuffers.release();
}

__global__ void massPropertiesKernel(
    const float3* __restrict__ positions,
    const float*  __restrict__ particleMasses,
    const int*    __restrict__ nodeDepth,
    const int*    __restrict__ nodeParticle,
    const int*    __restrict__ nodeChildren,
    const int*    __restrict__ nodeHasChildren,
    float*        __restrict__ nodeMass,
    float3*       __restrict__ nodeCenterOfMass,
    int nodeCount,
    int currentDepth)
{
    const int node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= nodeCount) return;
    if (nodeDepth[node] != currentDepth) return;

    if (!nodeHasChildren[node])
    {
        const int p = nodeParticle[node];
        if (p >= 0)
        {
            nodeMass[node] = particleMasses[p];
            nodeCenterOfMass[node] = positions[p];
        }
        else
        {
            nodeMass[node] = 0.f;
            nodeCenterOfMass[node] = make_float3(0.f, 0.f, 0.f);
        }
        return;
    }

    float totalMass = 0.f;
    float3 com = make_float3(0.f, 0.f, 0.f);
    const int base = node * 8;
    #pragma unroll
    for (int c = 0; c < 8; ++c)
    {
        const int child = nodeChildren[base + c];
        if (child < 0) continue;
        const float m = nodeMass[child];
        totalMass += m;
        com.x += nodeCenterOfMass[child].x * m;
        com.y += nodeCenterOfMass[child].y * m;
        com.z += nodeCenterOfMass[child].z * m;
    }
    nodeMass[node] = totalMass;
    if (totalMass > 0.f)
    {
        const float inv = 1.f / totalMass;
        com.x *= inv;
        com.y *= inv;
        com.z *= inv;
    }
    nodeCenterOfMass[node] = com;
}

__global__ void barnesHutForceKernel(
    const float3* __restrict__ positions,
    const float3* __restrict__ nodeCenterOfMass,
    const float*  __restrict__ nodeMass,
    const float*  __restrict__ nodeHalfSize,
    const int*    __restrict__ nodeChildren,
    const int*    __restrict__ nodeParticle,
    const int*    __restrict__ nodeHasChildren,
    float3*       __restrict__ accelerations,
    int   nParticles,
    int   nodeCount,
    float G,
    float theta2,
    float softening2)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nParticles) return;
    const float3 pi = positions[i];
    float3 acc = make_float3(0.f, 0.f, 0.f);

    int stack[STACK_CAP];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0)
    {
        const int node = stack[--sp];
        if (node < 0 || node >= nodeCount) continue;
        const float mass = nodeMass[node];
        if (mass <= 0.f) continue;
        const float3 com = nodeCenterOfMass[node];
        const float dx = com.x - pi.x;
        const float dy = com.y - pi.y;
        const float dz = com.z - pi.z;
        const float r2 = dx*dx + dy*dy + dz*dz + softening2;
        if (!nodeHasChildren[node])
        {
            const int other = nodeParticle[node];
            if (other == EMPTY || other == i) continue;
            const float inv_r  = rsqrtf(r2);
            const float inv_r3 = inv_r * inv_r * inv_r;
            const float factor = G * mass * inv_r3;
            acc.x += dx * factor;
            acc.y += dy * factor;
            acc.z += dz * factor;
            continue;
        }
        const float size = nodeHalfSize[node] * 2.f;
        const float size2 = size * size;
        if (size2 < theta2 * r2)
        {
            const float inv_r  = rsqrtf(r2);
            const float inv_r3 = inv_r * inv_r * inv_r;
            const float factor = G * mass * inv_r3;
            acc.x += dx * factor;
            acc.y += dy * factor;
            acc.z += dz * factor;
        }
        else
        {
            const int base = node * 8;
            #pragma unroll
            for (int c = 0; c < 8; ++c)
            {
                const int child = nodeChildren[base + c];
                if (child != EMPTY && sp < STACK_CAP)
                    stack[sp++] = child;
            }
        }
    }
    accelerations[i] = acc;
}


std::vector<Vector3> computeBarnesHutStepCuda(
    const std::vector<Vector3>& positions,
    const std::vector<float>&   particleMasses,
    const std::vector<float>&   nodeHalfSize,
    const std::vector<int>&     nodeDepth,
    const std::vector<int>&     nodeChildren,
    const std::vector<int>&     nodeParticle,
    const std::vector<bool>&    nodeHasChildren,
    int   nodeCount,
    int   maxDepth,
    float gravitationalConstant,
    float theta,
    float softening,
    int   blockSize)
{
    const int n = static_cast<int>(positions.size());
    if (n == 0) return {};
    if (nodeCount <= 0)
        throw std::runtime_error("BarnesHut CUDA: empty tree");
    if (static_cast<int>(particleMasses.size()) != n)
        throw std::runtime_error("positions / masses size mismatch");
    int block = blockSize;
    if (block < 32) block = 32;
    if (block > 512) block = 512;
    const float softening2 = softening * softening;
    const float theta2 = theta * theta;
    std::vector<float3> hPos(n);
    for (int i = 0; i < n; ++i)
        hPos[i] = make_float3(positions[i].x, positions[i].y, positions[i].z);
    std::vector<int> hHasChildren(nodeCount);
    for (int i = 0; i < nodeCount; ++i)
        hHasChildren[i] = nodeHasChildren[i] ? 1 : 0;
    gBuffers.ensure(n, nodeCount);
    const size_t childrenBytes = static_cast<size_t>(nodeCount) * 8 * sizeof(int);

    CUDA_CHECK(cudaMemcpy(gBuffers.dPos,      hPos.data(),            n * sizeof(float3), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gBuffers.dPartMass, particleMasses.data(),  n * sizeof(float),  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gBuffers.dHalfSize, nodeHalfSize.data(),    nodeCount * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gBuffers.dDepth,    nodeDepth.data(),       nodeCount * sizeof(int),   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gBuffers.dChildren, nodeChildren.data(),    childrenBytes,             cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gBuffers.dParticle, nodeParticle.data(),    nodeCount * sizeof(int),   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gBuffers.dHasChildren, hHasChildren.data(), nodeCount * sizeof(int),   cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(gBuffers.dNodeMass, 0, nodeCount * sizeof(float)));
    CUDA_CHECK(cudaMemset(gBuffers.dCom,      0, nodeCount * sizeof(float3)));
    const int nodeBlocks = (nodeCount + block - 1) / block;
    for (int depth = maxDepth; depth >= 0; --depth)
    {
        massPropertiesKernel<<<nodeBlocks, block>>>(
            gBuffers.dPos, gBuffers.dPartMass, gBuffers.dDepth, gBuffers.dParticle,
            gBuffers.dChildren, gBuffers.dHasChildren,
            gBuffers.dNodeMass, gBuffers.dCom,
            nodeCount, depth);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    const int particleBlocks = (n + block - 1) / block;
    barnesHutForceKernel<<<particleBlocks, block>>>(
        gBuffers.dPos, gBuffers.dCom, gBuffers.dNodeMass, gBuffers.dHalfSize,
        gBuffers.dChildren, gBuffers.dParticle, gBuffers.dHasChildren,
        gBuffers.dAcc, n, nodeCount, gravitationalConstant, theta2, softening2);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float3> hAcc(n);
    CUDA_CHECK(cudaMemcpy(hAcc.data(), gBuffers.dAcc, n * sizeof(float3), cudaMemcpyDeviceToHost));
    std::vector<Vector3> result(n);
    for (int i = 0; i < n; ++i)
        result[i] = Vector3(hAcc[i].x, hAcc[i].y, hAcc[i].z);
    return result;
}
}