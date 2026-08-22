#include "BarnesHutTree.h"
#include "particles/ParticleGroup.h"

#ifdef GALAXY_ENABLE_CUDA
#include "cuda/BarnesHutKernel.cuh"
#endif

#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <functional>


BarnesHutTree::BarnesHutTree(int maxParticles): maxParticles(maxParticles), maxNodes(8 * maxParticles), nodeCount(0) {
    nodeCenter.resize(maxNodes);
    nodeHalfSize.resize(maxNodes);
    nodeParent.resize(maxNodes, EMPTY);
    nodeDepth.resize(maxNodes, 0);
    nodeParticle.resize(maxNodes, EMPTY);
    nodeChildren.resize(maxNodes * 8, EMPTY);
    nodeHasChildren.resize(maxNodes, false);
    nodeMass.resize(maxNodes, 0.0f);
    nodeCenterOfMass.resize(maxNodes);
    clear();
}

void BarnesHutTree::clear()
{
    nodeCount = 0;
    std::fill(nodeParent.begin(), nodeParent.end(), EMPTY);
    std::fill(nodeDepth.begin(), nodeDepth.end(), 0);
    std::fill(nodeParticle.begin(), nodeParticle.end(), EMPTY);
    std::fill(nodeChildren.begin(), nodeChildren.end(), EMPTY);
    std::fill(nodeMass.begin(), nodeMass.end(), 0.0f);
    std::fill(nodeHalfSize.begin(), nodeHalfSize.end(), 0.0f);
    std::fill(nodeCenter.begin(), nodeCenter.end(), Vector3());
    std::fill(nodeCenterOfMass.begin(), nodeCenterOfMass.end(), Vector3());
    std::fill(nodeHasChildren.begin(), nodeHasChildren.end(), false);
}

void BarnesHutTree::build(const ParticleGroup& particles)
{
    clear();
    computeRootBounds(particles);
    for (int particle = 0; particle < particles.nbParticles; particle++) {insertParticle(particles, particle);}
    computeMassProperties(particles);
}

int BarnesHutTree::createNode()
{
    if (nodeCount >= maxNodes) {throw std::runtime_error("BarnesHutTree capacity exceeded.");}
    int node = nodeCount++;
    nodeCenter[node] = Vector3();
    nodeCenterOfMass[node] = Vector3();
    nodeHalfSize[node] = 0.0f;
    nodeMass[node] = 0.0f;
    nodeParent[node] = EMPTY;
    nodeDepth[node] = 0;
    nodeParticle[node] = EMPTY;
    for (int child = 0; child < 8; child++) {nodeChildren[node * 8 + child] = EMPTY;}
    return node;
}



int BarnesHutTree::childIndex(int node, const Vector3& position) const
{
    int child = 0;
    if (position.x >= nodeCenter[node].x) {child |= 1;}
    if (position.y >= nodeCenter[node].y) {child |= 2;}
    if (position.z >= nodeCenter[node].z) {child |= 4;}
    return child;
}



Vector3 BarnesHutTree::childCenter(int node, int child) const
{
    float offset = nodeHalfSize[node] * 0.5f;
    Vector3 center = nodeCenter[node];
    center.x += (child & 1) ? offset : -offset;
    center.y += (child & 2) ? offset : -offset;
    center.z += (child & 4) ? offset : -offset;
    return center;
}

void BarnesHutTree::computeRootBounds(const ParticleGroup& particles)
{
    if (particles.nbParticles == 0)
    {
        createNode();
        return;
    }
    float minX = std::numeric_limits<float>::max();
    float minY = std::numeric_limits<float>::max();
    float minZ = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::lowest();
    float maxY = std::numeric_limits<float>::lowest();
    float maxZ = std::numeric_limits<float>::lowest();

    for (int i = 0; i < particles.nbParticles; i++)
    {
        const Vector3& position = particles.positions[i];
        minX = std::min(minX, position.x);
        minY = std::min(minY, position.y);
        minZ = std::min(minZ, position.z);
        maxX = std::max(maxX, position.x);
        maxY = std::max(maxY, position.y);
        maxZ = std::max(maxZ, position.z);
    }
    float sizeX = maxX - minX;
    float sizeY = maxY - minY;
    float sizeZ = maxZ - minZ;
    float cubeSize = std::max({sizeX, sizeY, sizeZ});
    if (cubeSize == 0.0f) {cubeSize = 1.0f;}
    int root = createNode();
    nodeCenter[root] = Vector3(0.5f * (minX + maxX), 0.5f * (minY + maxY), 0.5f * (minZ + maxZ));
    nodeHalfSize[root] = 0.5f * cubeSize;
}

void BarnesHutTree::insertParticle(const ParticleGroup& particles, int particleIndex)
{
    int node = ROOT;
    while (true)
    {
        if (!nodeHasChildren[node] && nodeParticle[node] == EMPTY)
        {
            nodeParticle[node] = particleIndex;
            return;
        }
        if (!nodeHasChildren[node])
        {
            int existingParticle = nodeParticle[node];
            nodeParticle[node] = EMPTY;
            nodeHasChildren[node] = true;
            for (int child = 0; child < 8; child++)
            {
                int childNode = createNode();
                nodeChildren[node * 8 + child] = childNode;
                nodeParent[childNode] = node;
                nodeDepth[childNode] = nodeDepth[node] + 1;
                nodeHalfSize[childNode] = nodeHalfSize[node] * 0.5f;
                nodeCenter[childNode] = childCenter(node, child);
            }
            insertParticle(particles, existingParticle);
            continue;
        }
        int child = childIndex(node, particles.positions[particleIndex]);
        node = nodeChildren[node * 8 + child];
    }
}

void BarnesHutTree::computeMassProperties(const ParticleGroup& particles)
{
    if (nodeCount == 0) {return;}
    int maximumDepth = 0;
    for (int node = 0; node < nodeCount; node++) {maximumDepth = std::max(maximumDepth, nodeDepth[node]);}

#ifdef GALAXY_ENABLE_CUDA
    if (particles.device == "gpu") {return;}
#endif

    for (int node = 0; node < nodeCount; node++)
    {
        if (!nodeHasChildren[node])
        {
            int particle = nodeParticle[node];
            if (particle != EMPTY)
            {
                nodeMass[node] = particles.masses[particle];
                nodeCenterOfMass[node] = particles.positions[particle];
            }
            else
            {
                nodeMass[node] = 0.0f;
                nodeCenterOfMass[node] = Vector3();
            }
        }
        else
        {
            nodeMass[node] = 0.0f;
            nodeCenterOfMass[node] = Vector3();
        }
    }
    for (int depth = maximumDepth; depth >= 0; depth--)
    {
        for (int node = 0; node < nodeCount; node++)
        {
            if (nodeDepth[node] != depth) {continue;}
            if (!nodeHasChildren[node]) {continue;}
            float totalMass = 0.0f;
            Vector3 centerOfMass;
            for (int child = 0; child < 8; child++)
            {
                int childNode = nodeChildren[node * 8 + child];
                if (childNode == EMPTY) {continue;}
                float childMass = nodeMass[childNode];
                totalMass += childMass;
                centerOfMass.x += nodeCenterOfMass[childNode].x * childMass;
                centerOfMass.y += nodeCenterOfMass[childNode].y * childMass;
                centerOfMass.z += nodeCenterOfMass[childNode].z * childMass;
            }
            nodeMass[node] = totalMass;
            if (totalMass > 0.0f)
            {
                centerOfMass.x /= totalMass;
                centerOfMass.y /= totalMass;
                centerOfMass.z /= totalMass;
            }
            nodeCenterOfMass[node] =  centerOfMass;
        }
    }
}

std::vector<Vector3> BarnesHutTree::computeAccelerations(const ParticleGroup& particles, float gravitationalConstant, float theta, float softening) const
{
#ifdef GALAXY_ENABLE_CUDA
    if (particles.device == "gpu") {
        int maximumDepth = 0;
        for (int node = 0; node < nodeCount; ++node)
            maximumDepth = std::max(maximumDepth, nodeDepth[node]);
        const int blockSize = 256;
        return galaxy_cuda::computeBarnesHutStepCuda(
                    particles.positions,
                    particles.masses,
                    nodeHalfSize,
                    nodeDepth,
                    nodeChildren,
                    nodeParticle,
                    nodeHasChildren,
                    nodeCount,
                    maximumDepth,
                    gravitationalConstant,
                    theta,
                    softening,
                    blockSize);
    }
#endif

    std::vector<Vector3> accelerations(particles.nbParticles, Vector3());
    float softeningSquared = softening * softening;
    for (int particle = 0; particle < particles.nbParticles; particle++)
    {
        Vector3 acceleration;
        std::vector<int> stack;
        stack.push_back(ROOT);
        while (!stack.empty())
        {
            int node = stack.back();
            stack.pop_back();
            if (node == EMPTY) {continue;}
            if (!nodeHasChildren[node])
            {
                int otherParticle = nodeParticle[node];
                if (otherParticle == EMPTY || otherParticle == particle) {continue;}
                Vector3 displacement = particles.positions[otherParticle] - particles.positions[particle];
                float distanceSquared = displacement.x * displacement.x + displacement.y * displacement.y + displacement.z * displacement.z + softeningSquared;
                float inverseDistance = 1.0f / std::sqrt(distanceSquared);
                float inverseDistanceCubed = inverseDistance * inverseDistance * inverseDistance;
                float factor = gravitationalConstant * particles.masses[otherParticle] * inverseDistanceCubed;
                acceleration.x += displacement.x * factor;
                acceleration.y += displacement.y * factor;
                acceleration.z += displacement.z * factor;
                continue;
            }
            Vector3 displacement = nodeCenterOfMass[node] - particles.positions[particle];
            float distanceSquared = displacement.x * displacement.x + displacement.y * displacement.y + displacement.z * displacement.z + softeningSquared;
            float distance = std::sqrt(distanceSquared);
            float size = nodeHalfSize[node] * 2.0f;
            if (size / distance < theta)
            {
                float inverseDistance =  1.0f / distance;
                float inverseDistanceCubed = inverseDistance * inverseDistance * inverseDistance;
                float factor = gravitationalConstant * nodeMass[node] * inverseDistanceCubed;
                acceleration.x += displacement.x * factor;
                acceleration.y += displacement.y * factor;
                acceleration.z += displacement.z * factor;
            }
            else
            {
                for (int child = 0; child < 8; child++)
                {
                    int childNode = nodeChildren[node * 8 + child];
                    if (childNode != EMPTY) {stack.push_back(childNode);}
                }
            }
        }


        accelerations[particle] =  acceleration;
    }
    return accelerations;
}

int BarnesHutTree::getMaxParticles() const {return maxParticles;}