#include "BarnesHutKernel.cuh"

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>

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
constexpr int BLOCK     = 256;


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
    float theta,
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
        if (!nodeHasChildren[node])
        {
            const int other = nodeParticle[node];
            if (other == EMPTY || other == i) continue;
            const float3 com = nodeCenterOfMass[node];
            const float dx = com.x - pi.x;
            const float dy = com.y - pi.y;
            const float dz = com.z - pi.z;
            const float r2 = dx*dx + dy*dy + dz*dz + softening2;
            const float inv_r  = rsqrtf(r2);
            const float inv_r3 = inv_r * inv_r * inv_r;
            const float factor = G * mass * inv_r3;
            acc.x += dx * factor;
            acc.y += dy * factor;
            acc.z += dz * factor;
            continue;
        }
        const float3 com = nodeCenterOfMass[node];
        const float dx = com.x - pi.x;
        const float dy = com.y - pi.y;
        const float dz = com.z - pi.z;
        const float r2 = dx*dx + dy*dy + dz*dz + softening2;
        const float dist = sqrtf(r2);
        const float size = nodeHalfSize[node] * 2.f;
        if (size / dist < theta)
        {
            const float inv_r  = 1.f / dist;
            const float inv_r3 = inv_r * inv_r * inv_r;
            const float factor = G * mass * inv_r3;
            acc.x += dx * factor;
            acc.y += dy * factor;
            acc.z += dz * factor;
        }
        else
        {
            const int base = node * 8;
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

std::vector<Vector3> computeBarnesHutAccelerationsCuda(
    const std::vector<Vector3>& positions,
    const std::vector<Vector3>& nodeCenterOfMass,
    const std::vector<float>&   nodeMass,
    const std::vector<float>&   nodeHalfSize,
    const std::vector<int>&     nodeChildren,
    const std::vector<int>&     nodeParticle,
    const std::vector<bool>&    nodeHasChildren,
    int   nodeCount,
    float gravitationalConstant,
    float theta,
    float softening)
{
    const int n = static_cast<int>(positions.size());
    if (n == 0) return {};
    if (nodeCount <= 0)
        throw std::runtime_error("BarnesHut CUDA: empty tree");
    const float softening2 = softening * softening;
    std::vector<float3> hPos(n);
    for (int i = 0; i < n; ++i)
        hPos[i] = make_float3(positions[i].x, positions[i].y, positions[i].z);
    std::vector<float3> hCom(nodeCount);
    for (int i = 0; i < nodeCount; ++i)
        hCom[i] = make_float3(nodeCenterOfMass[i].x, nodeCenterOfMass[i].y, nodeCenterOfMass[i].z);
    std::vector<int> hHasChildren(nodeCount);
    for (int i = 0; i < nodeCount; ++i)
        hHasChildren[i] = nodeHasChildren[i] ? 1 : 0;
    float3* dPos          = nullptr;
    float3* dCom          = nullptr;
    float*  dMass         = nullptr;
    float*  dHalfSize     = nullptr;
    int*    dChildren     = nullptr;
    int*    dParticle     = nullptr;
    int*    dHasChildren  = nullptr;
    float3* dAcc          = nullptr;

    const size_t childrenBytes = static_cast<size_t>(nodeCount) * 8 * sizeof(int);
    CUDA_CHECK(cudaMalloc(&dPos,         n * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&dCom,         nodeCount * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&dMass,        nodeCount * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dHalfSize,    nodeCount * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dChildren,    childrenBytes));
    CUDA_CHECK(cudaMalloc(&dParticle,    nodeCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dHasChildren, nodeCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dAcc,         n * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(dPos,         hPos.data(),              n * sizeof(float3), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dCom,         hCom.data(),              nodeCount * sizeof(float3), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dMass,        nodeMass.data(),          nodeCount * sizeof(float),  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dHalfSize,    nodeHalfSize.data(),      nodeCount * sizeof(float),  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dChildren,    nodeChildren.data(),      childrenBytes,              cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dParticle,    nodeParticle.data(),      nodeCount * sizeof(int),    cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dHasChildren, hHasChildren.data(),      nodeCount * sizeof(int),    cudaMemcpyHostToDevice));

    const int blocks = (n + BLOCK - 1) / BLOCK;
    barnesHutForceKernel<<<blocks, BLOCK>>>(dPos, dCom, dMass, dHalfSize, dChildren, dParticle, dHasChildren, dAcc, n, nodeCount, gravitationalConstant, theta, softening2);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    std::vector<float3> hAcc(n);
    CUDA_CHECK(cudaMemcpy(hAcc.data(), dAcc, n * sizeof(float3), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(dPos));
    CUDA_CHECK(cudaFree(dCom));
    CUDA_CHECK(cudaFree(dMass));
    CUDA_CHECK(cudaFree(dHalfSize));
    CUDA_CHECK(cudaFree(dChildren));
    CUDA_CHECK(cudaFree(dParticle));
    CUDA_CHECK(cudaFree(dHasChildren));
    CUDA_CHECK(cudaFree(dAcc));
    std::vector<Vector3> result(n);
    for (int i = 0; i < n; ++i)
        result[i] = Vector3(hAcc[i].x, hAcc[i].y, hAcc[i].z);
    return result;
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
        // Leaf
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

void computeMassPropertiesCuda(
    const std::vector<Vector3>& positions,
    const std::vector<float>&   masses,
    const std::vector<int>&     nodeDepth,
    const std::vector<int>&     nodeParticle,
    const std::vector<int>&     nodeChildren,
    const std::vector<bool>&    nodeHasChildren,
    int nodeCount,
    int maxDepth,
    std::vector<float>&   nodeMass,
    std::vector<Vector3>& nodeCenterOfMass)
{
    if (nodeCount <= 0) return;

    const int nParticles = static_cast<int>(positions.size());
    std::vector<float3> hPos(nParticles);
    for (int i = 0; i < nParticles; ++i)
        hPos[i] = make_float3(positions[i].x, positions[i].y, positions[i].z);
    std::vector<int> hHasChildren(nodeCount);
    for (int i = 0; i < nodeCount; ++i)
        hHasChildren[i] = nodeHasChildren[i] ? 1 : 0;
    float3* dPos          = nullptr;
    float*  dPartMass     = nullptr;
    int*    dDepth        = nullptr;
    int*    dParticle     = nullptr;
    int*    dChildren     = nullptr;
    int*    dHasChildren  = nullptr;
    float*  dNodeMass     = nullptr;
    float3* dCom          = nullptr;
    const size_t childrenBytes = static_cast<size_t>(nodeCount) * 8 * sizeof(int);

    CUDA_CHECK(cudaMalloc(&dPos,         nParticles * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&dPartMass,    nParticles * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dDepth,       nodeCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dParticle,    nodeCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dChildren,    childrenBytes));
    CUDA_CHECK(cudaMalloc(&dHasChildren, nodeCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dNodeMass,    nodeCount * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dCom,         nodeCount * sizeof(float3)));
    CUDA_CHECK(cudaMemcpy(dPos,         hPos.data(),           nParticles * sizeof(float3), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dPartMass,    masses.data(),         nParticles * sizeof(float),  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dDepth,       nodeDepth.data(),      nodeCount * sizeof(int),     cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dParticle,    nodeParticle.data(),   nodeCount * sizeof(int),     cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dChildren,    nodeChildren.data(),   childrenBytes,               cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dHasChildren, hHasChildren.data(),   nodeCount * sizeof(int),     cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(dNodeMass, 0, nodeCount * sizeof(float)));
    CUDA_CHECK(cudaMemset(dCom,      0, nodeCount * sizeof(float3)));

    const int blocks = (nodeCount + BLOCK - 1) / BLOCK;
    for (int depth = maxDepth; depth >= 0; --depth)
    {
        massPropertiesKernel<<<blocks, BLOCK>>>(dPos, dPartMass, dDepth, dParticle, dChildren, dHasChildren, dNodeMass, dCom, nodeCount, depth);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    std::vector<float>  hMass(nodeCount);
    std::vector<float3> hCom(nodeCount);
    CUDA_CHECK(cudaMemcpy(hMass.data(), dNodeMass, nodeCount * sizeof(float),  cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hCom.data(),  dCom,      nodeCount * sizeof(float3), cudaMemcpyDeviceToHost));
    nodeMass.resize(nodeCount);
    nodeCenterOfMass.resize(nodeCount);
    for (int i = 0; i < nodeCount; ++i)
    {
        nodeMass[i] = hMass[i];
        nodeCenterOfMass[i] = Vector3(hCom[i].x, hCom[i].y, hCom[i].z);
    }
    CUDA_CHECK(cudaFree(dPos));
    CUDA_CHECK(cudaFree(dPartMass));
    CUDA_CHECK(cudaFree(dDepth));
    CUDA_CHECK(cudaFree(dParticle));
    CUDA_CHECK(cudaFree(dChildren));
    CUDA_CHECK(cudaFree(dHasChildren));
    CUDA_CHECK(cudaFree(dNodeMass));
    CUDA_CHECK(cudaFree(dCom));
}
    
}