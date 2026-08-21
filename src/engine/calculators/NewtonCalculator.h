#pragma once

#include "calculators/Calculator.h"


class NewtonCalculator : public Calculator
{
public:

    explicit NewtonCalculator(float gravitationalConstant = 1.0f, int tileSize = 16);

    std::vector<Vector3> computeAccelerations(const ParticleGroup& particles) override;

private:
    int tileSize;
};