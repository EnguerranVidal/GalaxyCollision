#pragma once

#include <vector>
#include "math/Vector3.h"

namespace galaxy_cuda
{
    std::vector<Vector3> computeNewtonAccelerationsCuda(
        const std::vector<Vector3>& positions,
        const std::vector<float>& masses,
        float gravitationalConstant,
        float softening,
        int tileSize);

    bool isCudaAvailable();
}