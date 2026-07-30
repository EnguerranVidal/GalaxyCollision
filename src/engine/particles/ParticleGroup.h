#pragma once

#include <vector>
#include <string>

#include "math/Vector3.h"


class ParticleGroup
{
public:
    explicit ParticleGroup(int nbParticles, const std::string& device = "cpu");

    void setPositions(const std::vector<Vector3>& positions);

    void setVelocities(const std::vector<Vector3>& velocities);

    void setMasses(const std::vector<float>& masses);

    void addParticle(const Vector3& position, const Vector3& velocity, float mass);

    [[nodiscard]] float kineticEnergy() const;

    void massCenter(Vector3& position, Vector3& velocity) const;

    [[nodiscard]] ParticleGroup groupToCpu() const;

    [[nodiscard]] ParticleGroup groupToGpu() const;

    [[nodiscard]] int getNbParticles() const;

    [[nodiscard]] const std::string& getDevice() const;

public:
    std::string device;
    int nbParticles = 0;
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;
    std::vector<Vector3> accelerations;
    std::vector<float> masses;
};