import numpy as np
from abc import ABC, abstractmethod

from src.core.engine.particles import ParticleGroup
from src.core.engine.kernels import massKernel, forceKernel

try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False
    cp = None


class Calculator(ABC):
    def __init__(self, gravitationalConstant: float = 1.0):
        self.gravitationalConstant = gravitationalConstant
        self.softening = 1e-3

    @abstractmethod
    def computeAccelerations(self, particles):
        pass


class NewtonCalculator(Calculator):
    def __init__(self, gravitationalConstant=1.0, blockSize=None):
        super().__init__(gravitationalConstant=gravitationalConstant)
        self.blockSize = blockSize

    def computeAccelerations(self, particles: ParticleGroup):
        xp = cp if particles.device == 'gpu' else np
        positions, masses = particles.positions, particles.masses
        nbParticles = particles.nbParticles
        accelerations = xp.empty_like(positions)
        blockSize = self._resolveBlockSize(nbParticles, particles.device)
        softeningSquared = self.softening * self.softening

        for startIndex in range(0, nbParticles, blockSize):
            endIndex = min(startIndex + blockSize, nbParticles)
            targetPositions = positions[startIndex:endIndex]
            displacements = positions[None, :, :] - targetPositions[:, None, :]
            distanceSquared = xp.sum(displacements * displacements, axis=2) + softeningSquared
            inverseDistanceCubed = distanceSquared ** -1.5
            weightedDisplacements = displacements * inverseDistanceCubed[:, :, None] * masses[None, :, None]
            accelerations[startIndex:endIndex] = self.gravitationalConstant * xp.sum(weightedDisplacements, axis=1)

        return accelerations

    def _resolveBlockSize(self, nbParticles, device):
        if self.blockSize is not None:
            return max(1, min(int(self.blockSize), nbParticles))
        if device == "gpu":
            return min(2048, nbParticles)
        return min(1024, nbParticles)


class BarnesHutCalculator(Calculator):
    def __init__(self, gravitationalConstant: float = 1.0, theta=0.5):
        super().__init__(gravitationalConstant=gravitationalConstant)
        self.theta = theta
        self.tree = None

    def computeAccelerations(self, particles):
        self._buildTree(particles)
        return self.tree.computeAccelerations(gravitationalConstant=self.gravitationalConstant, theta=self.theta, softening=self.softening)

    def _buildTree(self, particles):
        self.tree = BarnesHutTree(particles)
        self.tree.computeMassCenters()


class BarnesHutTree:
    ROOT = 0
    EMPTY = -1

    def __init__(self, particles):
        self.particles = particles
        xp = cp if particles.device == "gpu" else np
        self.xp = xp
        self.maxNodes = 4 * particles.nbParticles
        self.nodeCount = 1
        self.nodeCenter = self.xp.zeros((self.maxNodes, 3), dtype=self.xp.float32)
        self.nodeHalfSize = self.xp.zeros(self.maxNodes, dtype=self.xp.float32)
        self.nodeChildren = self.xp.full((self.maxNodes, 8), -1, dtype=self.xp.int32)
        self.nodeParent = self.xp.full(self.maxNodes, self.EMPTY, dtype=self.xp.int32)
        self.nodeDepth = self.xp.zeros(self.maxNodes, dtype=self.xp.int16)
        self.nodeChildCount = self.xp.zeros(self.maxNodes, dtype=self.xp.uint8)
        self.nodeParticle = self.xp.full(self.maxNodes, -1, dtype=self.xp.int32)
        self.nodeMass = self.xp.zeros(self.maxNodes, dtype=self.xp.float32)
        self.nodeCenterOfMass = self.xp.zeros((self.maxNodes, 3), dtype=self.xp.float32)
        self.nodeIsLeaf = self.xp.ones(self.maxNodes, dtype=bool)
        self.build()

    @property
    def positions(self):
        return self.particles.positions

    @property
    def masses(self):
        return self.particles.masses

    @property
    def nbParticles(self):
        return self.particles.nbParticles

    def build(self):
        self._initializeRoot()
        for particleIndex in range(self.nbParticles):
            self._insertParticle(particleIndex)

    def _initializeRoot(self):
        minimum = self.xp.min(self.positions, axis=0)
        maximum = self.xp.max(self.positions, axis=0)
        self.nodeCenter[self.ROOT] = 0.5 * (minimum + maximum)
        half = float(self.xp.max(maximum - minimum)) * 0.5
        if half <= 0:
            half = 1e-6
        self.nodeHalfSize[self.ROOT] = half

    @staticmethod
    def _octant(position, center):
        octant = 0
        if position[0] >= center[0]:
            octant |= 1
        if position[1] >= center[1]:
            octant |= 2
        if position[2] >= center[2]:
            octant |= 4
        return octant

    def _allocateNode(self):
        if self.nodeCount >= self.maxNodes:
            raise RuntimeError("Barnes-Hut tree capacity exceeded.")
        node = self.nodeCount
        self.nodeCount += 1
        return node

    def _childCenter(self, parent, octant):
        half = self.nodeHalfSize[parent] * 0.5
        center = self.nodeCenter[parent].copy()
        center[0] += half if (octant & 1) else -half
        center[1] += half if (octant & 2) else -half
        center[2] += half if (octant & 4) else -half
        return center, half

    def _createChild(self, parent, octant):
        child = self.nodeChildren[parent, octant]
        if child != self.EMPTY:
            return child
        child = self._allocateNode()
        self.nodeChildren[parent, octant] = child
        self.nodeParent[child] = parent
        self.nodeDepth[child] = self.nodeDepth[parent] + 1
        self.nodeChildCount[parent] += 1
        center, half = self._childCenter(parent, octant)
        self.nodeCenter[child] = center
        self.nodeHalfSize[child] = half
        return child

    def _subdivide(self, node):
        self.nodeIsLeaf[node] = False

    def _insertParticle(self, particleIndex):
        current = self.ROOT
        while True:
            if self.nodeIsLeaf[current]:
                storedParticle = self.nodeParticle[current]
                if storedParticle == self.EMPTY:
                    self.nodeParticle[current] = particleIndex
                    return
                self._subdivide(current)
                self.nodeParticle[current] = self.EMPTY
                oldOctant = self._octant(self.positions[storedParticle], self.nodeCenter[current])
                child = self._createChild(current, oldOctant)
                self.nodeParticle[child] = storedParticle
            octant = self._octant(self.positions[particleIndex], self.nodeCenter[current])
            current = self._createChild(current, octant)

    def computeMassCenters(self):
        if self.particles.device == "gpu":
            self._computeMassCentersGpu()
        else:
            self._computeMassCentersCpu()

    def _computeMassCentersCpu(self):
        self.nodeMass[:self.nodeCount] = 0
        self.nodeCenterOfMass[:self.nodeCount] = 0
        for node in range(self.nodeCount):
            if self.nodeIsLeaf[node]:
                particle = self.nodeParticle[node]
                if particle != self.EMPTY:
                    self.nodeMass[node] = self.masses[particle]
                    self.nodeCenterOfMass[node] = self.positions[particle]
        for node in range(self.nodeCount - 1, -1, -1):
            if not self.nodeIsLeaf[node]:
                totalMass = 0.0
                center = self.xp.zeros(3, dtype=self.xp.float32)
                for child in self.nodeChildren[node]:
                    if child != self.EMPTY:
                        mass = self.nodeMass[child]
                        totalMass += mass
                        center += mass * self.nodeCenterOfMass[child]
                if totalMass > 0:
                    self.nodeMass[node] = totalMass
                    self.nodeCenterOfMass[node] = center / totalMass

    def _computeMassCentersGpu(self):
        self.nodeMass.fill(0)
        self.nodeCenterOfMass.fill(0)
        maxDepth = int(cp.max(self.nodeDepth[:self.nodeCount]).get())
        for depth in range(maxDepth, -1, -1):
            massKernel(..., np.int32(depth), np.int32(self.nodeCount), block=(256,), grid=((self.nodeCount + 255) // 256,))

    def computeAccelerations(self, gravitationalConstant, theta, softening):
        if self.particles.device == "gpu":
            return self._computeAccelerationsGpu(gravitationalConstant, theta, softening)
        return self._computeAccelerationsCpu(gravitationalConstant, theta, softening)

    def _computeAccelerationsCpu(self, gravitationalConstant=1.0, theta=0.5, softening=1e-3):
        accelerations = self.xp.zeros_like(self.positions)
        softeningSquared = softening * softening
        for particleIndex in range(self.nbParticles):
            position = self.positions[particleIndex]
            acceleration = self.xp.zeros(3, dtype=self.xp.float32)
            stack = [self.ROOT]
            while stack:
                node = stack.pop()
                if self.nodeMass[node] == 0:
                    continue
                if self.nodeIsLeaf[node] and self.nodeParticle[node] == particleIndex:
                    continue
                displacement = self.nodeCenterOfMass[node] - position
                distanceSquared = self.xp.dot(displacement, displacement) + softeningSquared
                distance = self.xp.sqrt(distanceSquared)
                size = self.nodeHalfSize[node] * 2
                if self.nodeIsLeaf[node] or size / distance < theta:
                    inverseDistanceCubed = distanceSquared ** -1.5
                    acceleration += gravitationalConstant * self.nodeMass[node] * displacement * inverseDistanceCubed
                else:
                    for child in self.nodeChildren[node]:
                        if child != self.EMPTY:
                            stack.append(int(child))
            accelerations[particleIndex] = acceleration
        return accelerations

    def _computeAccelerationsGpu(self, gravitationalConstant=1.0, theta=0.5, softening=1e-3):
        accelerations = cp.zeros_like(self.positions)
        threads = 256
        blocks = (self.nbParticles + threads - 1) // threads
        self._forceKernel(self.positions, self.nodeCenterOfMass, self.nodeMass, self.nodeHalfSize, self.nodeChildren, self.nodeParticle, self.nodeIsLeaf, np.int32(self.nodeCount), np.float32(gravitationalConstant), np.float32(theta), np.float32(softening), accelerations, block=(threads,), grid=(blocks,))
        return accelerations
