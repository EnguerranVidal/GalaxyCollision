#pragma once

#include "math/Vector3.h"
#include "particles/ParticleGroup.h"
#include "parameters/SimulatorParameters.h"


class Distribution
{

public:

    virtual ~Distribution() = default;
    virtual ParticleGroup generate(const SimulatorParameters& parameters) = 0;

protected:
    static ParticleGroup createParticleGroup(const std::vector<Vector3>& positions, const std::vector<Vector3>& velocities, const std::vector<float>& masses, const std::string& device);

};