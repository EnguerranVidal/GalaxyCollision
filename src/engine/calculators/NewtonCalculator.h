#pragma once

#include "calculators/Calculator.h"


class NewtonCalculator : public Calculator
{
public:

    explicit NewtonCalculator(float gravitationalConstant = 1.0f);

    std::vector<Vector3> computeAccelerations(const ParticleGroup& particles) override;
};