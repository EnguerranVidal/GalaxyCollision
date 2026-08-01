#pragma once

#include "Distribution.h"


class BasicDistribution : public Distribution
{
public:
    ParticleGroup generate(const SimulatorParameters& parameters) override;
};