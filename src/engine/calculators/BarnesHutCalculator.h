#pragma once

#include <memory>

#include "calculators/Calculator.h"
#include "calculators/BarnesHutTree.h"


class BarnesHutCalculator : public Calculator
{
public:

    BarnesHutCalculator(float gravitationalConstant = 1.0f, float theta = 0.5f);

    std::vector<Vector3> computeAccelerations(const ParticleGroup& particles) override;

private:

    void buildTree(const ParticleGroup& particles );

    float theta;
    std::unique_ptr<BarnesHutTree> tree;
};