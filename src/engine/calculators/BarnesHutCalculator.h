#pragma once

#include <memory>

#include "calculators/Calculator.h"
#include "calculators/BarnesHutTree.h"


class BarnesHutCalculator : public Calculator
{
public:

    BarnesHutCalculator(float gravitationalConstant = 1.0f, float theta = 0.5f, int blockSize = 256);

    std::vector<Vector3> computeAccelerations(const ParticleGroup& particles) override;

private:

    void buildTree(const ParticleGroup& particles );

    float theta;
    int blockSize;
    std::unique_ptr<BarnesHutTree> tree;
};