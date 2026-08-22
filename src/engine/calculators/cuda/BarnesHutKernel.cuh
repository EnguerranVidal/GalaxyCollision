#pragma once

#include <vector>
#include "math/Vector3.h"

namespace galaxy_cuda
{
    std::vector<Vector3> computeBarnesHutStepCuda(
        const std::vector<Vector3>& positions,
        const std::vector<float>&   particleMasses,
        const std::vector<float>&   nodeHalfSize,
        const std::vector<int>&     nodeDepth,
        const std::vector<int>&     nodeChildren,   // node*8 + child, EMPTY=-1
        const std::vector<int>&     nodeParticle,
        const std::vector<bool>&    nodeHasChildren,
        int   nodeCount,
        int   maxDepth,
        float gravitationalConstant,
        float theta,
        float softening,
        int   blockSize = 256);

    void releaseBarnesHutCudaBuffers();

}