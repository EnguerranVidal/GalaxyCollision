import numpy as np

try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False
    cp = None


class ParticleGroup:
        def __init__(self, nbParticles, device='cpu'):
            xp = cp if device == 'gpu' else np
            self.nbParticles = nbParticles
            self.device = device
            self.positions = xp.zeros((self.nbParticles, 3), dtype=xp.float32)
            self.velocities = xp.zeros((self.nbParticles, 3), dtype=xp.float32)
            self.accelerations = xp.zeros((self.nbParticles, 3), dtype=xp.float32)
            self.masses = xp.ones(self.nbParticles, dtype=xp.float32)

        def setPositions(self, positions):
            self.positions = cp.array(positions, dtype=cp.float32) if self.device == 'gpu' else np.array(positions, dtype=np.float32)

        def setVelocities(self, velocities):
            self.velocities = cp.array(velocities, dtype=cp.float32) if self.device == 'gpu' else np.array(velocities, dtype=np.float32)

        def setMasses(self, masses):
            self.masses = cp.array(masses, dtype=cp.float32) if self.device == 'gpu' else np.array(masses, dtype=np.float32)

        def addParticle(self, position, velocity, mass):
            xp = cp if self.device == 'gpu' else np
            singlePosition = xp.array([position], dtype=xp.float32).reshape(1, 3)
            singleVelocity = xp.array([velocity], dtype=xp.float32).reshape(1, 3)
            singleMass = xp.array([mass], dtype=xp.float32)
            self.positions = xp.concatenate([self.positions, singlePosition])
            self.velocities = xp.concatenate([self.velocities, singleVelocity])
            self.masses = xp.concatenate([self.masses, singleMass])
            self.accelerations = xp.concatenate([self.accelerations, xp.zeros((1, 3), dtype=xp.float32)])
            self.nbParticles += 1

        def kineticEnergy(self):
            xp = cp if self.device == 'gpu' else np
            squaredVelocity = xp.sum(self.velocities ** 2, axis=1)
            return 0.5 * xp.sum(self.masses * squaredVelocity)

        def massCenter(self):
            xp = cp if self.device == 'gpu' else np
            totalMass = xp.sum(self.masses)
            massCenterPosition = xp.sum(self.masses[:, None] * self.positions, axis=0) / totalMass
            massCenterVelocity = xp.sum(self.masses[:, None] * self.velocities, axis=0) / totalMass
            return massCenterPosition, massCenterVelocity

        def groupToCpu(self):
            if self.device == 'cpu':
                return self
            particleGroup = ParticleGroup(self.nbParticles, device='cpu')
            particleGroup.positions = self.positions.get()
            particleGroup.velocities = self.velocities.get()
            particleGroup.accelerations = self.accelerations.get()
            particleGroup.masses = self.masses.get()
            return particleGroup

        def groupToGpu(self):
            if self.device == 'gpu':
                return self
            particleGroup = ParticleGroup(self.nbParticles, device='gpu')
            particleGroup.positions = cp.array(self.positions, dtype=cp.float32)
            particleGroup.velocities = cp.array(self.velocities, dtype=cp.float32)
            particleGroup.accelerations = cp.array(self.accelerations, dtype=cp.float32)
            particleGroup.masses = cp.array(self.masses, dtype=cp.float32)
            return particleGroup
