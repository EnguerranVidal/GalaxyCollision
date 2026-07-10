import numpy as np


class Initializer:
    def __init__(self, gravitationalConstant=1.0):
        self.gravitationalConstant = gravitationalConstant

    @staticmethod
    def basicDistribution(numParticles: int, massRange=(0.1, 1.0), positionScale=10.0, velocityScale=1.0, is3D=True):
        dim = 3 if is3D else 2
        positions = (np.random.rand(numParticles, dim) * 2 - 1) * positionScale
        masses = np.random.uniform(massRange[0], massRange[1], numParticles)
        velocities = (np.random.randn(numParticles, dim) * velocityScale)
        return positions, velocities, masses

    def galaxyDistribution(self, numParticles: int, totalMass=1000.0, radius=10.0, height=2.0, is3D=True):
        if is3D:
            r = np.random.exponential(scale=radius / 2, size=numParticles)
            theta = np.random.uniform(0, 2 * np.pi, numParticles)
            z = np.random.normal(0, height / 4, numParticles)
            x, y = r * np.cos(theta), r * np.sin(theta)
            positions = np.stack([x, y, z], axis=1)
        else:
            r = np.random.exponential(scale=radius / 2, size=numParticles)
            theta = np.random.uniform(0, 2 * np.pi, numParticles)
            x, y = r * np.cos(theta), r * np.sin(theta)
            positions = np.stack([x, y], axis=1)
        masses = np.full(numParticles, totalMass / numParticles, dtype=np.float64)
        if is3D:
            distance = np.sqrt(np.sum(positions[:, :2] ** 2, axis=1) + 1e-6)
            circularVelocity = np.sqrt(self.gravitationalConstant * totalMass / distance)
            xVelocity = -positions[:, 1] / distance * circularVelocity
            yVelocity = positions[:, 0] / distance * circularVelocity
            zVelocities = np.random.normal(0, 0.1, numParticles)
            velocities = np.stack([xVelocity, yVelocity, zVelocities], axis=1)
        else:
            distance = np.sqrt(np.sum(positions ** 2, axis=1) + 1e-6)
            circularVelocity = np.sqrt(self.gravitationalConstant * totalMass / distance)
            xVelocity = -positions[:, 1] / distance * circularVelocity
            yVelocity = positions[:, 0] / distance * circularVelocity
            velocities = np.stack([xVelocity, yVelocity], axis=1)
        velocities += np.random.normal(0, 0.2, velocities.shape)
        return positions, velocities, masses