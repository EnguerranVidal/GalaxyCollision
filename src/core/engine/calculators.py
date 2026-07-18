import cupy as cp
import numpy as np
from abc import ABC, abstractmethod

from src.core.engine.particles import ParticleGroup


class Calculator(ABC):
    def __init__(self, gravitationalConstant: float = 1.0):
        self.gravitationalConstant = gravitationalConstant
        self.softening = 1e-3

    @abstractmethod
    def computeAccelerations(self, particles):
        pass


class NewtonCalculator(Calculator):
    def __init__(self, gravitationalConstant=1.0):
        super().__init__(gravitationalConstant=gravitationalConstant)

    def computeAccelerations(self, particles: ParticleGroup):
        positions, velocities, masses = particles.positions.copy(), particles.velocities.copy(), particles.masses.copy()
        xp = cp if particles.device == 'gpu' else np
        massMatrix = masses.reshape((1, -1, 1)) * masses.reshape((-1, 1, 1))
        displacements = positions.reshape((1, -1, 3)) - positions.reshape((-1, 1, 3))
        distances = xp.linalg.norm(displacements, axis=2)
        forces = self.gravitationalConstant * displacements * massMatrix / xp.expand_dims(distances + self.softening, 2) ** 3
        return forces.sum(axis=1) / masses.reshape(-1, 1)
