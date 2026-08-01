#include "Distribution.h"


ParticleGroup Distribution::createParticleGroup(const std::vector<Vector3>& positions, const std::vector<Vector3>& velocities, const std::vector<float>& masses, const std::string& device)
{
    ParticleGroup particleGroup(static_cast<int>(positions.size()), device);
    particleGroup.setPositions(positions);
    particleGroup.setVelocities(velocities);
    particleGroup.setMasses(masses);
    return particleGroup;
}