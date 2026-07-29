import cupy as cp
from cupyx import jit

@jit.rawkernel()
def massKernel(nodeDepth, nodeIsLeaf, nodeParticle, nodeChildren, nodeMass, nodeCenterOfMass, particleMass, particlePositions, currentDepth, nodeCount):
    node = jit.blockIdx.x * jit.blockDim.x + jit.threadIdx.x
    if node >= nodeCount:
        return
    if nodeDepth[node] != currentDepth:
        return
    if nodeIsLeaf[node]:
        particle = nodeParticle[node]
        if particle >= 0:
            nodeMass[node] = particleMass[particle]
            nodeCenterOfMass[node, 0] = particlePositions[particle, 0]
            nodeCenterOfMass[node, 1] = particlePositions[particle, 1]
            nodeCenterOfMass[node, 2] = particlePositions[particle, 2]
        return
    totalMass = 0.0
    centerX = 0.0
    centerY = 0.0
    centerZ = 0.0
    count = 0
    for childIndex in range(8):
        child = nodeChildren[node, childIndex]
        if child == -1:
            continue
        count += 1
        if count == nodeCount[node]:
            break
        mass = nodeMass[child]
        totalMass += mass
        centerX += mass * nodeCenterOfMass[child, 0]
        centerY += mass * nodeCenterOfMass[child, 1]
        centerZ += mass * nodeCenterOfMass[child, 2]
    nodeMass[node] = totalMass
    if totalMass > 0:
        inverseMass = 1.0 / totalMass
        nodeCenterOfMass[node, 0] = centerX * inverseMass
        nodeCenterOfMass[node, 1] = centerY * inverseMass
        nodeCenterOfMass[node, 2] = centerZ * inverseMass


@jit.rawkernel()
def forceKernel(positions, nodeCenterOfMass, nodeMass, nodeHalfSize, nodeChildren, nodeParticle, nodeIsLeaf, nodeCount, gravitationalConstant, theta, softeningSquared, accelerations):
    particle = jit.blockIdx.x * jit.blockDim.x + jit.threadIdx.x
    if particle >= positions.shape[0]:
        return
    px = positions[particle,0]
    py = positions[particle,1]
    pz = positions[particle,2]
    ax = 0.0
    ay = 0.0
    az = 0.0
    stack = jit.local_array((64,), dtype=cp.int32)
    stackSize = 0
    stack[stackSize] = 0
    stackSize += 1
    while stackSize > 0:
        stackSize -= 1
        node = stack[stackSize]
        mass = nodeMass[node]
        if mass <= 0:
            continue
        if nodeIsLeaf[node]:
            if nodeParticle[node] == particle:
                continue
        dx = nodeCenterOfMass[node,0] - px
        dy = nodeCenterOfMass[node,1] - py
        dz = nodeCenterOfMass[node,2] - pz
        distanceSquared = dx * dx + dy * dy + dz * dz + softeningSquared
        size = nodeHalfSize[node] * 2.0
        accept = nodeIsLeaf[node] or size * size < theta * theta * distanceSquared
        if accept:
            inverseDistance = 1.0 / jit.rsqrt(distanceSquared)
            inverseDistanceCubed = inverseDistance * inverseDistance * inverseDistance
            factor = gravitationalConstant * mass * inverseDistanceCubed
            ax += factor * dx
            ay += factor * dy
            az += factor * dz
        else:
            for childIndex in range(8):
                child = nodeChildren[node, childIndex]
                if child != -1:
                    if stackSize < 64:
                        stack[stackSize] = child
                        stackSize += 1
    accelerations[particle,0] = ax
    accelerations[particle,1] = ay
    accelerations[particle,2] = az
