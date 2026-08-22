#include "GalaxyDistribution.h"

#include <random>
#include <cmath>
#include <numbers>
#include <algorithm>

constexpr float PI = std::numbers::pi_v<float>;


ParticleGroup GalaxyDistribution::generate(const SimulatorParameters& parameters)
{
    const auto& p = parameters.galaxyDistributionParameters;
    const int nbParticles = parameters.nbParticles;
    std::mt19937 generator(parameters.seed);
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;
    std::vector<float> masses;
    positions.reserve(nbParticles);
    velocities.reserve(nbParticles);
    masses.reserve(nbParticles);
    generateDisk(p, nbParticles, generator, positions, velocities, masses);
    generateBulge(p, nbParticles, generator, positions, velocities, masses);
    generateHalo(p, nbParticles,generator, positions, velocities, masses);
    return createParticleGroup(positions, velocities, masses, parameters.device);
}

void GalaxyDistribution::generateDisk(const GalaxyDistributionParameters& p, int nbParticles, std::mt19937& generator, std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<float>& masses)
{
    int nbDisk = static_cast<int>(static_cast<float>(nbParticles) * (1.0f - p.bulgeFraction - p.haloFraction));
    std::exponential_distribution<float> exponential(1.0f);
    std::uniform_real_distribution<float> thetaDistribution(0.0f, 2.0f * PI);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    float diskMass = p.totalMass * (1.0f - p.bulgeFraction - p.haloFraction);
    for(int i = 0; i < nbDisk; i++)
    {
        float radius = p.radius * exponential(generator);
        float theta = thetaDistribution(generator);
        float height = normal(generator) * p.height / 4.0f;
        positions.emplace_back(radius * cos(theta), radius * sin(theta), height);
        float enclosedMass = diskMass * (1.0f - (1.0f + radius / p.radius) * exp(-radius / p.radius));
        float circularVelocity = sqrt(std::max(enclosedMass / std::max(radius,0.01f), 0.0f));
        float sigma = p.velocityDispersion * circularVelocity;
        float tangential = normal(generator) * sigma;
        float radial = normal(generator) * sigma;
        velocities.emplace_back(
            -tangential * sin(theta) + radial * cos(theta),
            tangential * cos(theta) + radial * sin(theta),
            normal(generator) * sigma * 0.5f
        );
        masses.emplace_back(diskMass / static_cast<float>(nbDisk));
    }
}

void GalaxyDistribution::generateBulge(const GalaxyDistributionParameters& p, int nbParticles, std::mt19937& generator, std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<float>& masses)
{
    int nbBulge = static_cast<int>(static_cast<float>(nbParticles) * p.bulgeFraction);
    std::uniform_real_distribution<float> uniform(0.01f, 1.0f);
    std::uniform_real_distribution<float> cosDistribution(-1.0f, 1.0f);
    std::uniform_real_distribution<float> phiDistribution(0.0f, 2.0f * PI);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    float bulgeMass = p.totalMass * p.bulgeFraction;
    for(int i=0;i<nbBulge;i++)
    {
        float u = uniform(generator);
        float radius = std::min( p.plummerRadius / sqrt(pow(u,-2.0f / 3.0f) - 1.0f), 10.0f * p.plummerRadius);
        float cosTheta = cosDistribution(generator);
        float phi = phiDistribution(generator);
        float sinTheta = sqrt(1.0f-cosTheta*cosTheta);
        positions.emplace_back(radius*sinTheta*cos(phi), radius*sinTheta*sin(phi), radius*cosTheta);
        float sigma = sqrt(std::max(bulgeMass / (6.0f * sqrt(radius * radius + p.plummerRadius * p.plummerRadius)), 0.0f));
        velocities.emplace_back(normal(generator) * sigma, normal(generator) * sigma, normal(generator) * sigma);
        masses.emplace_back(bulgeMass / static_cast<float>(nbBulge));
    }
}

void GalaxyDistribution::generateHalo(const GalaxyDistributionParameters& p, int nbParticles, std::mt19937& generator, std::vector<Vector3>& positions, std::vector<Vector3>& velocities, std::vector<float>& masses)
{
    int nbHalo = nbParticles - static_cast<int>(static_cast<float>(nbParticles) * (1.0f - p.bulgeFraction - p.haloFraction)) - static_cast<int>(static_cast<float>(nbParticles) * p.bulgeFraction);
    std::exponential_distribution<float> exponential(1.0f);
    std::uniform_real_distribution<float> cosDistribution(-1.0f, 1.0f);
    std::uniform_real_distribution<float> phiDistribution(0.0f, 2.0f * PI);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    float haloMass = p.totalMass * p.haloFraction;
    for(int i=0;i<nbHalo;i++)
    {
        float radius = std::min(p.haloRadius * exponential(generator), 20.0f * p.haloRadius);
        float cosTheta = cosDistribution(generator);
        float phi = phiDistribution(generator);
        float sinTheta = sqrt(1.0f-cosTheta*cosTheta);
        positions.emplace_back(radius * sinTheta * cos(phi), radius * sinTheta * sin(phi), radius * cosTheta);
        float sigma = sqrt(haloMass / (radius + p.haloRadius));
        velocities.emplace_back(normal(generator) * sigma, normal(generator) * sigma, normal(generator) * sigma);
        masses.emplace_back(haloMass / static_cast<float>(nbHalo));
    }
}