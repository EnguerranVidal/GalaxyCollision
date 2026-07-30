#pragma once

#include "integrators/Integrator.h"
#include "particles/ParticleGroup.h"
#include "calculators/Calculator.h"


class EulerExplicitIntegrator : public Integrator
{
public:
    explicit EulerExplicitIntegrator(float timeStep);

    void step(ParticleGroup& particles, Calculator& calculator) override;
};