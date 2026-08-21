#include "NewtonCalculator.h"
#include "particles/ParticleGroup.h"

#ifdef GALAXY_ENABLE_CUDA
#include "cuda/NewtonKernel.cuh"
#endif

#include <cmath>
#include <stdexcept>


NewtonCalculator::NewtonCalculator(float gravitationalConstant, int tileSize): Calculator(gravitationalConstant), tileSize(tileSize) {}

std::vector<Vector3> NewtonCalculator::computeAccelerations(const ParticleGroup& particles)
{
#ifdef GALAXY_ENABLE_CUDA
    if (particles.device == "gpu")
    {
        if (!galaxy_cuda::isCudaAvailable())
            throw std::runtime_error("GPU device requested but no CUDA device is available");
        return galaxy_cuda::computeNewtonAccelerationsCuda(particles.positions, particles.masses, gravitationalConstant, softening, tileSize);
    }
#endif

    std::vector<Vector3> accelerations(particles.nbParticles, Vector3());

    const float softeningSquared = softening * softening;
    for (int i = 0; i < particles.nbParticles; i++)
    {
        Vector3 acceleration;
        const Vector3& positionI = particles.positions[i];
        for (int j = 0; j < particles.nbParticles; j++)
        {
            if (i == j) {continue;}
            const Vector3& positionJ = particles.positions[j];
            float dx = positionJ.x - positionI.x;
            float dy = positionJ.y - positionI.y;
            float dz = positionJ.z - positionI.z;
            float distanceSquared = dx * dx + dy * dy + dz * dz + softeningSquared;
            float inverseDistance = 1.0f / std::sqrt(distanceSquared);
            float inverseDistanceCubed = inverseDistance * inverseDistance * inverseDistance;
            float factor = gravitationalConstant * particles.masses[j] * inverseDistanceCubed;
            acceleration.x += dx * factor;
            acceleration.y += dy * factor;
            acceleration.z += dz * factor;
        }
        accelerations[i] = acceleration;
    }
    return accelerations;
}