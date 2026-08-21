#include "BarnesHutCalculator.h"
#include "particles/ParticleGroup.h"


BarnesHutCalculator::BarnesHutCalculator(float gravitationalConstant, float theta, int blockSize): Calculator(gravitationalConstant), theta(theta), blockSize(blockSize) {}

std::vector<Vector3> BarnesHutCalculator::computeAccelerations(const ParticleGroup& particles)
{
    buildTree(particles);
    return tree->computeAccelerations(particles, gravitationalConstant, theta, softening);
}

void BarnesHutCalculator::buildTree(const ParticleGroup& particles)
{
    if (!tree || tree->getMaxParticles() < particles.nbParticles) {tree = std::make_unique<BarnesHutTree>(particles.nbParticles);}
    tree->build(particles);
}