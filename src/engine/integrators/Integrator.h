#pragma once

#include "particles/ParticleGroup.h"
#include "calculators/Calculator.h"


class Integrator
{
public:

    explicit Integrator(float timeStep);

    virtual ~Integrator() = default;

    virtual void step(ParticleGroup& particles, Calculator& calculator) = 0;

    [[nodiscard]] float getTimeStep() const;

protected:
    float timeStep;
};