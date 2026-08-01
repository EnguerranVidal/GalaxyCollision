#include "BasicDistribution.h"

#include <random>


ParticleGroup BasicDistribution::generate(const SimulatorParameters& parameters)
{
    const auto& distributionParameters = parameters.basicDistributionParameters;
    std::mt19937 generator(parameters.seed);
    std::uniform_real_distribution<float> positionDistribution(-1.0f, 1.0f);
    std::uniform_real_distribution<float> massDistribution(distributionParameters.massMinimum, distributionParameters.massMaximum);
    std::normal_distribution<float> velocityDistribution(0.0f, 1.0f);
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;
    std::vector<float> masses;
    positions.reserve(distributionParameters.nbParticles);
    velocities.reserve(distributionParameters.nbParticles);
    masses.reserve(distributionParameters.nbParticles);
    for(int i = 0; i < distributionParameters.nbParticles; i++)
    {
        positions.emplace_back(
            positionDistribution(generator) * distributionParameters.positionScale,
            positionDistribution(generator) * distributionParameters.positionScale,
            positionDistribution(generator) * distributionParameters.positionScale
        );
        velocities.emplace_back(
            velocityDistribution(generator) * distributionParameters.velocityScale,
            velocityDistribution(generator) * distributionParameters.velocityScale,
            velocityDistribution(generator) * distributionParameters.velocityScale
        );
        masses.emplace_back(massDistribution(generator));
    }
    return createParticleGroup(positions, velocities, masses, parameters.device);
}