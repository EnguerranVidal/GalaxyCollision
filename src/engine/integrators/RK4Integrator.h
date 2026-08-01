#pragma once

#include "integrators/Integrator.h"
#include "particles/ParticleGroup.h"
#include "calculators/Calculator.h"


class RK4Integrator : public Integrator
{
public:
    explicit RK4Integrator(float timeStep);

    void step(ParticleGroup& particles, Calculator& calculator) override;
};