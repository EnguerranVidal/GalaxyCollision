#pragma once

#include <vector>

#include "math/Vector3.h"
#include "particles/ParticleGroup.h"


class Calculator
{
public:

    explicit Calculator(float gravitationalConstant = 1.0f);

    virtual ~Calculator() = default;

    virtual std::vector<Vector3> computeAccelerations(const ParticleGroup& particles) = 0;

    [[nodiscard]] float getGravitationalConstant() const;

    [[nodiscard]] float getSoftening() const;

protected:
    float gravitationalConstant;
    float softening;
};