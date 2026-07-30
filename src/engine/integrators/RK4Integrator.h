#pragma once

#include "integrators/Integrator.h"
#include "particles/ParticleGroup.h"
#include "calculators/Calculator.h"


class RK4 : public Integrator
{
public:
    explicit RK4(float timeStep);

    void step(ParticleGroup& particles, Calculator& calculator) override;
};