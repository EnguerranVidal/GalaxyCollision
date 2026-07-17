import numpy as np

from src.core.engine.parameters import BasicDistributionParameters, GalaxyDistributionParameters
from src.core.engine.particles import ParticleGroup


class Distribution:
    def __init__(self, gravitationalConstant=1.0, device='cpu'):
        self.device = device
        self.nbParticles = 0
        self.gravitationalConstant = gravitationalConstant
        self.positionOffset = np.zeros(3, dtype=np.float32)
        self.velocityOffset = np.zeros(3, dtype=np.float32)

    def setPosition(self, position):
        self.positionOffset = position

    def setVelocity(self, velocity):
        self.velocityOffset = velocity

    def generate(self, parameters: BasicDistributionParameters = None, seed=None):
        randomGenerator = np.random.RandomState(seed)
        positions = (randomGenerator.rand(parameters.nbParticles, 3) * 2 - 1) * parameters.positionScale
        masses = randomGenerator.uniform(parameters.massMinimum, parameters.massMaximum, parameters.nbParticles)
        velocities = randomGenerator.randn(parameters.nbParticles, 3) * parameters.velocityScale
        return self.createParticleGroup(positions, velocities, masses)

    def createParticleGroup(self, positions, velocities, masses):
        positions += self.positionOffset
        velocities += self.velocityOffset
        particleGroup = ParticleGroup(nbParticles=len(positions), device=self.device)
        particleGroup.setPositions(positions)
        particleGroup.setVelocities(velocities)
        particleGroup.setMasses(masses)
        return particleGroup

    def __len__(self):
        return self.nbParticles


class GalaxyDistribution(Distribution):
    def generate(self, parameters: GalaxyDistributionParameters, seed=None):
        randomGenerator = np.random.RandomState(seed)

        # EXPONENTIAL GALAXY DISK DISTRIBUTION
        nbDisk = int(parameters.nbParticles * (1 - parameters.bulgeFraction - parameters.haloFraction))
        diskRadii = parameters.radius * randomGenerator.exponential(1.0, nbDisk)
        diskTheta = randomGenerator.uniform(0, 2 * np.pi, nbDisk)
        diskHeight = randomGenerator.normal(0, parameters.height / 4, nbDisk)
        diskPositions = np.column_stack([diskRadii * np.cos(diskTheta), diskRadii * np.sin(diskTheta), diskHeight])
        enlacedMasses = parameters.totalMass * (1 - parameters.bulgeFraction - parameters.haloFraction) * (1 - (1 + diskRadii / parameters.radius) * np.exp(-diskRadii / parameters.radius))
        circularVelocities = np.sqrt(np.maximum(enlacedMasses / np.maximum(diskRadii, 0.01), 0))
        diskSigma = parameters.velocityDispersion * circularVelocities
        tangentialVelocities, radialVelocities = randomGenerator.normal(0, diskSigma) , randomGenerator.normal(0, diskSigma)
        diskVelocities = np.column_stack([-tangentialVelocities * np.sin(diskTheta) + radialVelocities * np.cos(diskTheta), tangentialVelocities * np.cos(diskTheta) + radialVelocities * np.sin(diskTheta), randomGenerator.normal(0, diskSigma * 0.5, nbDisk)])
        diskMasses = np.full(nbDisk, parameters.totalMass * (1 - parameters.bulgeFraction - parameters.haloFraction) / nbDisk)

        # PLUMMER SPHERICAL BULGE DISTRIBUTION
        nbBulge = int(parameters.nbParticles * parameters.bulgeFraction)
        uniformDistribution = randomGenerator.uniform(0.01, 1.0, nbBulge)
        bulgeRadii = np.minimum(parameters.plummerRadius / np.sqrt(uniformDistribution ** (-2.0 / 3.0) - 1), 10 * parameters.plummerRadius)
        bulgeCosTheta = randomGenerator.uniform(-1, 1, nbBulge)
        bulgePhi = randomGenerator.uniform(0, 2 * np.pi, nbBulge)
        bulgeSinTheta = np.sqrt(1 - bulgeCosTheta**2)
        bulgePositions = np.column_stack([bulgeRadii * bulgeSinTheta * np.cos(bulgePhi), bulgeRadii * bulgeSinTheta * np.sin(bulgePhi), bulgeRadii * bulgeCosTheta])
        totalBulgeMass = parameters.totalMass * parameters.bulgeFraction
        bulgeSigma = np.sqrt(np.maximum(totalBulgeMass / (6 * np.sqrt(bulgeRadii ** 2 + parameters.plummerRadius ** 2)), 0))
        bulgeVelocities = np.column_stack([randomGenerator.normal(0, bulgeSigma), randomGenerator.normal(0, bulgeSigma), randomGenerator.normal(0, bulgeSigma)])
        bulgeMasses = np.full(nbBulge, totalBulgeMass / nbBulge)

        # DARK MATTER HALO DISTRIBUTION
        nbHalo = parameters.nbParticles - nbDisk - nbBulge
        haloRadii = np.minimum(parameters.haloRadius * randomGenerator.exponential(1.0, nbHalo), 20 * parameters.haloRadius)
        haloCosTheta = randomGenerator.uniform(-1, 1, nbHalo)
        haloPhi = randomGenerator.uniform(0, 2 * np.pi, nbHalo)
        haloSinTheta = np.sqrt(1 - haloCosTheta ** 2)
        haloPositions = np.column_stack([haloRadii * haloSinTheta * np.cos(haloPhi), haloRadii * haloSinTheta * np.sin(haloPhi), haloRadii * haloCosTheta])
        totalHaloMass = parameters.totalMass * parameters.bulgeFraction
        haloSigma = np.sqrt(totalHaloMass / (haloRadii + parameters.haloRadius))
        haloVelocities = np.column_stack([randomGenerator.normal(0, haloSigma), randomGenerator.normal(0, haloSigma), randomGenerator.normal(0, haloSigma)])
        haloMasses = np.full(nbHalo, totalHaloMass / nbHalo)

        # CONCATENATE POSITION, VELOCITIES, MASSES
        positions = np.concatenate([diskPositions, bulgePositions, haloPositions])
        masses = np.concatenate([diskMasses, bulgeMasses, haloMasses])
        velocities = np.concatenate([diskVelocities, bulgeVelocities, haloVelocities])
        return self.createParticleGroup(positions, velocities, masses)
