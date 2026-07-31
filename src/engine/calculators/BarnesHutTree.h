#pragma once

#include <vector>

#include "math/Vector3.h"
#include "particles/ParticleGroup.h"


class BarnesHutTree
{
public:

    static constexpr int ROOT = 0;
    static constexpr int EMPTY = -1;

    explicit BarnesHutTree(int maxParticles);

    void build(const ParticleGroup& particles);

    void clear();


private:

    void computeRootBounds(const ParticleGroup& particles);


    void insertParticle(const ParticleGroup& particles, int particleIndex);

    void computeMassProperties(const ParticleGroup& particles);

    int createNode();

    [[nodiscard]] int childIndex(int node, const Vector3& position) const;

    [[nodiscard]] Vector3 childCenter(int node,int child) const;

    [[nodiscard]] std::vector<Vector3> computeAccelerations(const ParticleGroup& particles, float gravitationalConstant, float theta, float softening) const;

    int maxParticles;
    int maxNodes;
    int nodeCount;
    std::vector<Vector3> nodeCenter;
    std::vector<float> nodeHalfSize;
    std::vector<int> nodeParent;
    std::vector<int> nodeDepth;
    std::vector<int> nodeParticle;
    std::vector<int> nodeChildren;
    std::vector<bool> nodeHasChildren;
    std::vector<float> nodeMass;
    std::vector<Vector3> nodeCenterOfMass;
};