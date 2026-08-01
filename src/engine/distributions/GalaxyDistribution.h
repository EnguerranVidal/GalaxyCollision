#pragma once

#include "Distribution.h"
#include <random>


class GalaxyDistribution : public Distribution
{
public:
    ParticleGroup generate(const SimulatorParameters& parameters) override;

private:
    static void generateDisk(const GalaxyDistributionParameters& parameters, std::mt19937& generator, std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<float>& masses);

    static void generateBulge(const GalaxyDistributionParameters& parameters, std::mt19937& generator, std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<float>& masses);

    static void generateHalo(const GalaxyDistributionParameters& parameters, std::mt19937& generator, std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<float>& masses);
};