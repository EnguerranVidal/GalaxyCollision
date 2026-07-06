import cupy as cp


class Initializer:
    def __init__(self, gravitationalConstant=1.0):
        self.gravitationalConstant = gravitationalConstant

    @staticmethod
    def basicDistribution(numParticles: int, massRange=(0.1, 1.0), positionScale=10.0, velocityScale=1.0, is3D=True):
        dim = 3 if is3D else 2
        positions = (cp.random.rand(numParticles, dim) * 2 - 1) * positionScale
        masses = cp.random.uniform(massRange[0], massRange[1], numParticles)
        velocities = (cp.random.randn(numParticles, dim) * velocityScale)
        return positions, velocities, masses

    def galaxyDistribution(self, numParticles: int, totalMass=1000.0, radius=10.0, height=2.0, is3D=True):
        if is3D:
            r = cp.random.exponential(scale=radius / 2, size=numParticles)
            theta = cp.random.uniform(0, 2 * cp.pi, numParticles)
            z = cp.random.normal(0, height / 4, numParticles)
            x, y = r * cp.cos(theta), r * cp.sin(theta)
            positions = cp.stack([x, y, z], axis=1)
        else:
            r = cp.random.exponential(scale=radius / 2, size=numParticles)
            theta = cp.random.uniform(0, 2 * cp.pi, numParticles)
            x, y = r * cp.cos(theta), r * cp.sin(theta)
            positions = cp.stack([x, y], axis=1)
        masses = cp.full(numParticles, totalMass / numParticles, dtype=cp.float64)
        if is3D:
            distance = cp.sqrt(cp.sum(positions[:, :2] ** 2, axis=1) + 1e-6)
            circularVelocity = cp.sqrt(self.gravitationalConstant * totalMass / distance)
            xVelocity = -positions[:, 1] / distance * circularVelocity
            yVelocity = positions[:, 0] / distance * circularVelocity
            zVelocities = cp.random.normal(0, 0.1, numParticles)
            velocities = cp.stack([xVelocity, yVelocity, zVelocities], axis=1)
        else:
            distance = cp.sqrt(cp.sum(positions ** 2, axis=1) + 1e-6)
            circularVelocity = cp.sqrt(self.gravitationalConstant * totalMass / distance)
            xVelocity = -positions[:, 1] / distance * circularVelocity
            yVelocity = positions[:, 0] / distance * circularVelocity
            velocities = cp.stack([xVelocity, yVelocity], axis=1)
        velocities += cp.random.normal(0, 0.2, velocities.shape)
        return positions, velocities, masses