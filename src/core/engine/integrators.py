import numpy as np

from src.core.engine.calculators import Calculator
from src.core.engine.particles import ParticleGroup

try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False
    cp = None


class Integrator:
    def __init__(self, timeStep: float):
        self.timeStep = timeStep

    def step(self, particles: ParticleGroup, calculator: Calculator):
        raise NotImplementedError("Subclasses must implement step()")


class EulerExplicit(Integrator):
    def step(self, particles: ParticleGroup, calculator: Calculator):
        positions, velocities = particles.positions.copy(), particles.velocities.copy()
        accelerations = calculator.computeAccelerations(particles)
        particles.velocities = velocities + accelerations * self.timeStep
        particles.positions = positions + velocities * self.timeStep
        particles.accelerations = accelerations
        return particles


class RK4(Integrator):
    def step(self, particles, calculator):
        positions, velocities = particles.positions.copy(), particles.velocities.copy()
        # K1 COEFFICIENT
        k1Accelerations = calculator.computeAccelerations(particles)
        k1Velocities = k1Accelerations * self.timeStep
        k1Positions = velocities * self.timeStep
        # K2 COEFFICIENT
        particles.positions = positions + 0.5 * k1Positions
        particles.velocities = velocities + 0.5 * k1Velocities
        k2Accelerations = calculator.computeAccelerations(particles)
        k2Velocities = k2Accelerations * self.timeStep
        k2Positions = (velocities + 0.5 * k1Velocities) * self.timeStep
        # K3 COEFFICIENT
        particles.positions = positions + 0.5 * k2Positions
        particles.velocities = velocities + 0.5 * k2Velocities
        k3Accelerations = calculator.computeAccelerations(particles)
        k3Velocities = k3Accelerations * self.timeStep
        k3Positions = (velocities + 0.5 * k2Velocities) * self.timeStep
        # K4 COEFFICIENT
        particles.positions = positions + k3Positions
        particles.velocities = velocities + k3Velocities
        k4Acceleration = calculator.computeAccelerations(particles)
        k4Velocities = k4Acceleration * self.timeStep
        k4Positions = (velocities + k3Velocities) * self.timeStep
        # WEIGHTED AVERAGES
        particles.positions = positions + (k1Positions + 2 * k2Positions + 2 * k3Positions + k4Positions) / 6
        particles.velocities = velocities + (k1Velocities + 2 * k2Velocities + 2 * k3Velocities + k4Velocities) / 6
        particles.accelerations = k4Acceleration
        return particles
