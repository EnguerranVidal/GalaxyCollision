class Integrator:
    def __init__(self, timeStep: float):
        self.timeStep = timeStep

    def step(self, positions, velocities, masses, calculator):
        raise NotImplementedError("Subclasses must implement step()")


class EulerExplicit(Integrator):
    def step(self, positions, velocities, masses, calculator):
        acc = calculator.computeAccelerations(positions, masses)
        velocities += acc * self.timeStep
        positions += velocities * self.timeStep
        return positions, velocities


class RK4(Integrator):
    def step(self, positions, velocities, masses, calculator):
        # K1 COEFFICIENT
        k1Acceleration = calculator.computeAccelerations(positions, masses)
        k1Position = velocities
        k1Velocities = k1Acceleration
        # K2 COEFFICIENT
        k2Positions = positions + 0.5 * self.timeStep * k1Position
        k2Velocities = velocities + 0.5 * self.timeStep * k1Velocities
        k2Acceleration = calculator.computeAccelerations(k2Positions, masses)
        k2Position = k2Velocities
        k2Velocities = k2Acceleration
        # K3 COEFFICIENT
        k3Positions = positions + 0.5 * self.timeStep * k2Position
        k3Velocities = velocities + 0.5 * self.timeStep * k2Velocities
        k3Acceleration = calculator.computeAccelerations(k3Positions, masses)
        k3Position = k3Velocities
        k3Velocities = k3Acceleration
        # K4 COEFFICIENT
        k4Positions = positions + self.timeStep * k3Position
        k4Velocities = velocities + self.timeStep * k3Velocities
        k4Acceleration = calculator.computeAccelerations(k4Positions, masses)
        k4Positions = k4Velocities
        k4Velocities = k4Acceleration
        # WEIGHTED AVERAGES
        newPositions = positions + (self.timeStep / 6) * (k1Position + 2 * k2Position + 2 * k3Position + k4Positions)
        newVelocities = velocities + (self.timeStep / 6) * (k1Velocities + 2 * k2Velocities + 2 * k3Velocities + k4Velocities)
        return newPositions, newVelocities
