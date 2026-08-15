#include "NewtonKernel.cuh"

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>


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

bool isCudaAvailable()
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0);
}


constexpr int TILE_SIZE = 16;

__global__ void newtonForceTiledKernel(
    const float3* __restrict__ positions,
    const float*  __restrict__ masses,
    float3*       __restrict__ accelerations,
    int   n,
    float G,
    float softening2)
{
    __shared__ float3 shPos[TILE_SIZE];
    __shared__ float  shMass[TILE_SIZE];

    const int i = blockIdx.x * TILE_SIZE + threadIdx.x;
    float3 acc = make_float3(0.f, 0.f, 0.f);
    float3 pi  = (i < n) ? positions[i] : make_float3(0.f, 0.f, 0.f);

    // Number of tiles needed to cover all particles
    const int numTiles = (n + TILE_SIZE - 1) / TILE_SIZE;

    for (int tile = 0; tile < numTiles; ++tile)
    {
        // Load one tile of source particles into shared memory
        const int j = tile * TILE_SIZE + threadIdx.x;
        if (j < n)
        {
            shPos[threadIdx.x]  = positions[j];
            shMass[threadIdx.x] = masses[j];
        }
        else
        {
            shPos[threadIdx.x]  = make_float3(0.f, 0.f, 0.f);
            shMass[threadIdx.x] = 0.f;
        }
        __syncthreads();
        if (i < n)
        {
            #pragma unroll
            for (int k = 0; k < TILE_SIZE; ++k)
            {
                const int jj = tile * TILE_SIZE + k;
                if (jj >= n || jj == i) continue;

                const float3 pj = shPos[k];
                const float  dx = pj.x - pi.x;
                const float  dy = pj.y - pi.y;
                const float  dz = pj.z - pi.z;
                const float  r2 = dx*dx + dy*dy + dz*dz + softening2;

                const float inv_r  = rsqrtf(r2);
                const float inv_r3 = inv_r * inv_r * inv_r;
                const float factor = G * shMass[k] * inv_r3;

                acc.x += dx * factor;
                acc.y += dy * factor;
                acc.z += dz * factor;
            }
        }
        __syncthreads();
    }

    if (i < n)
        accelerations[i] = acc;
}


std::vector<Vector3> computeNewtonAccelerationsCuda(
    const std::vector<Vector3>& positions,
    const std::vector<float>& masses,
    float gravitationalConstant,
    float softening)
{
    const int n = static_cast<int>(positions.size());
    if (n == 0)
        return {};

    if (static_cast<int>(masses.size()) != n)
        throw std::runtime_error("positions and masses size mismatch");

    const float softening2 = softening * softening;

    std::vector<float3> hPos(n);
    for (int i = 0; i < n; ++i)
        hPos[i] = make_float3(positions[i].x, positions[i].y, positions[i].z);

    float3* dPos  = nullptr;
    float*  dMass = nullptr;
    float3* dAcc  = nullptr;

    CUDA_CHECK(cudaMalloc(&dPos,  n * sizeof(float3)));
    CUDA_CHECK(cudaMalloc(&dMass, n * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dAcc,  n * sizeof(float3)));

    CUDA_CHECK(cudaMemcpy(dPos,  hPos.data(),   n * sizeof(float3), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dMass, masses.data(), n * sizeof(float),  cudaMemcpyHostToDevice));

    const int blocks = (n + TILE_SIZE - 1) / TILE_SIZE;
    newtonForceTiledKernel<<<blocks, TILE_SIZE>>>(
        dPos, dMass, dAcc, n, gravitationalConstant, softening2);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<float3> hAcc(n);
    CUDA_CHECK(cudaMemcpy(hAcc.data(), dAcc, n * sizeof(float3), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(dPos));
    CUDA_CHECK(cudaFree(dMass));
    CUDA_CHECK(cudaFree(dAcc));

    std::vector<Vector3> result(n);
    for (int i = 0; i < n; ++i)
        result[i] = Vector3(hAcc[i].x, hAcc[i].y, hAcc[i].z);

    return result;
}

}