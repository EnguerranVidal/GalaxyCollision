#include "NewtonCalculator.h"

#include <cmath>

#include "particles/ParticleGroup.h"


NewtonCalculator::NewtonCalculator(float gravitationalConstant): Calculator(gravitationalConstant) {}

std::vector<Vector3> NewtonCalculator::computeAccelerations(const ParticleGroup& particles)
{
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