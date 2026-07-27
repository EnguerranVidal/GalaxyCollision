import numpy as np
from abc import ABC, abstractmethod

from src.core.engine.particles import ParticleGroup

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
