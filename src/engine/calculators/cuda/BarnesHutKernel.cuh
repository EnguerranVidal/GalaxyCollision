#pragma once

#include <vector>
#include "math/Vector3.h"

namespace galaxy_cuda
{
    std::vector<Vector3> computeBarnesHutAccelerationsCuda(
        const std::vector<Vector3>& positions,
        const std::vector<Vector3>& nodeCenterOfMass,
        const std::vector<float>&   nodeMass,
        const std::vector<float>&   nodeHalfSize,
        const std::vector<int>&     nodeChildren,
        const std::vector<int>&     nodeParticle,
        const std::vector<bool>&    nodeHasChildren,
        int   nodeCount,
        float gravitationalConstant,
        float theta,
        float softening);

    void computeMassPropertiesCuda(
        const std::vector<Vector3>& positions,
        const std::vector<float>&   masses,
        const std::vector<int>&     nodeDepth,
        const std::vector<int>&     nodeParticle,
        const std::vector<int>&     nodeChildren,
        const std::vector<bool>&    nodeHasChildren,
        int nodeCount,
        int maxDepth,
        std::vector<float>&   nodeMass,
        std::vector<Vector3>& nodeCenterOfMass);
}