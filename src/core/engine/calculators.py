import cupy as cp
from abc import ABC, abstractmethod


class Calculator(ABC):
    def __init__(self, gravitationalConstant: float = 1.0, is3D: bool = True):
        self.gravitationalConstant = gravitationalConstant
        self.is3D = is3D
        self.softening = 1e-3

    @abstractmethod
    def computeAccelerations(self, positions, masses):
        pass

    def _computePairwiseDistance(self, posI, posJ):
        dx = posJ - posI
        distSq = cp.sum(dx * dx, axis=-1) + self.softening ** 2
        return dx, cp.sqrt(distSq)


class NewtonCalculator(Calculator):
    def __init__(self, gravitationalConstant=1.0, is3D=True):
        super().__init__(gravitationalConstant=gravitationalConstant, is3D=is3D)

    def computeAccelerations(self, positions, masses):
        n = len(positions)
        dimension = positions.shape[1]
        accelerations = cp.zeros((n, dimension), dtype=cp.float64)
        for i in range(n):
            xDistance = positions[i + 1:] - positions[i]
            distSquared = cp.sum(xDistance ** 2, axis=1) + self.softening ** 2
            distance = cp.sqrt(distSquared)
            forceMagnitude = masses[i + 1:] / distSquared
            accelerations[i] += cp.sum(forceMagnitude[:, None] * xDistance / distance[:, None], axis=0)
            accelerations[i + 1:] -= (forceMagnitude[:, None] * xDistance / distance[:, None]) * (masses[i] / masses[i + 1:])[:, None]
        return accelerations * self.gravitationalConstant


class Node:
    def __init__(self, center, size, is3D=False):
        self.center = center
        self.size = size
        self.is3D = is3D
        self.mass = 0.0
        self.massCenter = cp.zeros(3 if is3D else 2, dtype=cp.float64)
        self.particles = []
        self.children = []
        self.isLeaf = True


class BarnesHutCalculator(Calculator):
    def __init__(self, gravitationalConstant=1.0, is3D=True, theta=0.5):
        super().__init__(gravitationalConstant=gravitationalConstant, is3D=is3D)
        self.theta = theta
        self.root = None

    def buildTree(self, positions, masses):
        n = len(positions)
        if n == 0:
            return
        minPosition = cp.min(positions, axis=0)
        maxPosition = cp.max(positions, axis=0)
        center = (minPosition + maxPosition) / 2.0
        size = cp.max(maxPosition - minPosition) * 1.1
        self.root = Node(center, size, self.is3D)
        for i in range(n):
            self._insert(self.root, positions[i], masses[i], i)

    def _insert(self, node, pos, mass, index):
        if node.isLeaf:
            if len(node.particles) == 0:
                node.particles.append(index)
                node.mass = mass
                node.centerOfMass = pos.copy()
                return
            node.isLeaf = False
            self._subdivide(node)
            oldPosition = node.centerOfMass
            oldMass = node.mass
            oldIndex = node.particles[0]
            node.particles = []
            self._insert(node, oldPosition, oldMass, oldIndex)
            self._insert(node, pos, mass, index)
        else:
            node.mass += mass
            node.centerOfMass = (node.centerOfMass * (node.mass - mass) + pos * mass) / node.mass
            childIndex = self._getChildIndex(node, pos)
            self._insert(node.children[childIndex], pos, mass, index)

    def _subdivide(self, node):
        half = node.size / 2.0
        offsets = cp.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]])[: (8 if self.is3D else 4)]
        for offset in offsets:
            childCenter = node.center + (offset - 0.5) * half
            child = Node(childCenter, half, self.is3D)
            node.children.append(child)

    def _getChildIndex(self, node, pos):
        idx = 0
        for dimension in range(3 if self.is3D else 2):
            if pos[dimension] > node.center[dimension]:
                idx |= (1 << dimension)
        return idx

    def computeAccelerations(self, positions, masses):
        self.buildTree(positions, masses)
        n = len(positions)
        accelerations = cp.zeros_like(positions)
        for i in range(n):
            accelerations[i] = self._calculateForce(self.root, positions[i], i)
        return accelerations * self.gravitationalConstant

    def _calculateForce(self, node, pos, particleIndex):
        if node.mass == 0:
            return cp.zeros_like(pos)
        xDistance = node.massCenter - pos
        distanceSquared = cp.dot(xDistance, xDistance) + self.softening ** 2
        distance = cp.sqrt(distanceSquared)
        if node.isLeaf or (node.size / distance < self.theta):
            if len(node.particles) == 1 and node.particles[0] == particleIndex:
                return cp.zeros_like(pos)
            forceMagnitude = node.mass / distanceSquared
            return forceMagnitude * xDistance / distance
        else:
            totalForce = cp.zeros_like(pos)
            for child in node.children:
                totalForce += self._calculateForce(child, pos, particleIndex)
            return totalForce
