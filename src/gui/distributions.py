import numpy as np

import engine
from src.gui.parameters import SimulatorParameters


def toParticleGroup(positions: np.ndarray, velocities: np.ndarray, masses: np.ndarray, device: str) -> engine.ParticleGroup:
    particleGroup = engine.ParticleGroup(positions.shape[0], device.lower())
    particleGroup.setPositions([engine.Vector3(float(x), float(y), float(z)) for x, y, z in positions])
    particleGroup.setVelocities([engine.Vector3(float(x), float(y), float(z)) for x, y, z in velocities])
    particleGroup.setMasses(masses.astype(np.float32).tolist())
    return particleGroup


def generateBasicDistribution(parameters: SimulatorParameters) -> engine.ParticleGroup:
    randomGenerator = np.random.default_rng(parameters.seed if parameters.seed is not None else None)
    n = parameters.nbParticles
    distributionParameters = parameters.basicDistributionParameters
    positions = (randomGenerator.uniform(-1.0, 1.0, size=(n, 3)) * distributionParameters.positionScale).astype(np.float32)
    velocities = (randomGenerator.normal(0.0, 1.0, size=(n, 3)) * distributionParameters.velocityScale).astype(np.float32)
    masses = randomGenerator.uniform(distributionParameters.massMinimum, distributionParameters.massMaximum, size=n).astype(np.float32)
    return toParticleGroup(positions, velocities, masses, parameters.device)


def generateGalaxyDistribution(parameters: SimulatorParameters) -> engine.ParticleGroup:
    randomGenerator = np.random.default_rng(parameters.seed if parameters.seed is not None else None)
    nbParticles = parameters.nbParticles
    distributionParameters = parameters.galaxyDistributionParameters
    nbDisk = int(nbParticles * (1.0 - distributionParameters.bulgeFraction - distributionParameters.haloFraction))
    nbBulge = int(nbParticles * distributionParameters.bulgeFraction)
    nbHalo = nbParticles - nbDisk - nbBulge
    positions = np.zeros((nbParticles, 3), dtype=np.float32)
    velocities = np.zeros((nbParticles, 3), dtype=np.float32)
    masses = np.zeros(nbParticles, dtype=np.float32)
    if nbDisk > 0:
        diskRadii = distributionParameters.radius * randomGenerator.exponential(1.0, size=nbDisk)
        diskTheta = randomGenerator.uniform(0.0, 2.0 * np.pi, size=nbDisk)
        diskHeight = randomGenerator.normal(0.0, distributionParameters.height / 4.0, size=nbDisk)
        positions[:nbDisk, 0] = diskRadii * np.cos(diskTheta)
        positions[:nbDisk, 1] = diskRadii * np.sin(diskTheta)
        positions[:nbDisk, 2] = diskHeight
        diskMass = distributionParameters.totalMass * (1.0 - distributionParameters.bulgeFraction - distributionParameters.haloFraction)
        masses[:nbDisk] = diskMass / nbDisk
        diskCircularVelocity = np.sqrt(np.maximum(diskMass / np.maximum(diskRadii, 0.01), 0.0))
        diskSigma = distributionParameters.velocityDispersion * diskCircularVelocity
        diskVt = diskCircularVelocity + randomGenerator.normal(0.0, diskSigma)
        diskVr = randomGenerator.normal(0.0, diskSigma)
        velocities[:nbDisk, 0] = -diskVt * np.sin(diskTheta) + diskVr * np.cos(diskTheta)
        velocities[:nbDisk, 1] = diskVt * np.cos(diskTheta) + diskVr * np.sin(diskTheta)
        velocities[:nbDisk, 2] = randomGenerator.normal(0.0, diskSigma * 0.5, size=nbDisk)
    if nbBulge > 0:
        uniformDistribution = randomGenerator.uniform(0.01, 1.0, nbBulge)
        bulgeRadii = np.minimum(distributionParameters.plummerRadius / np.sqrt(uniformDistribution ** (-2.0 / 3.0) - 1), 10 * distributionParameters.plummerRadius)
        bulgeCosTheta = randomGenerator.uniform(-1, 1, nbBulge)
        bulgePhi = randomGenerator.uniform(0, 2 * np.pi, nbBulge)
        bulgeSinTheta = np.sqrt(1 - bulgeCosTheta ** 2)
        positions[nbDisk: nbDisk + nbBulge, 0] = bulgeRadii * bulgeSinTheta * np.cos(bulgePhi)
        positions[nbDisk: nbDisk + nbBulge, 1] = bulgeRadii * bulgeSinTheta * np.sin(bulgePhi)
        positions[nbDisk: nbDisk + nbBulge, 2] = bulgeRadii * bulgeCosTheta
        totalBulgeMass = distributionParameters.totalMass * distributionParameters.bulgeFraction
        masses[nbDisk: nbDisk + nbBulge] = totalBulgeMass / nbBulge
        bulgeSigma = np.sqrt(np.maximum(totalBulgeMass / (6 * np.sqrt(bulgeRadii ** 2 + distributionParameters.plummerRadius ** 2)), 0))
        velocities[nbDisk: nbDisk + nbBulge, 0] = randomGenerator.normal(0, bulgeSigma)
        velocities[nbDisk: nbDisk + nbBulge, 1] = randomGenerator.normal(0, bulgeSigma)
        velocities[nbDisk: nbDisk + nbBulge, 2] = randomGenerator.normal(0, bulgeSigma)
    if nbHalo> 0:
        haloRadii = np.minimum(distributionParameters.haloRadius * randomGenerator.exponential(1.0, nbHalo), 20 * distributionParameters.haloRadius)
        haloCosTheta = randomGenerator.uniform(-1, 1, nbHalo)
        haloPhi = randomGenerator.uniform(0, 2 * np.pi, nbHalo)
        haloSinTheta = np.sqrt(1 - haloCosTheta ** 2)
        positions[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo, 0] = haloRadii * haloSinTheta * np.cos(haloPhi)
        positions[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo, 1] = haloRadii * haloSinTheta * np.sin(haloPhi)
        positions[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo, 2] = haloRadii * haloCosTheta
        totalHaloMass = distributionParameters.totalMass * distributionParameters.bulgeFraction
        masses[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo] = totalBulgeMass / nbBulge
        haloSigma = np.sqrt(totalHaloMass / (haloRadii + distributionParameters.haloRadius))
        velocities[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo, 0] = randomGenerator.normal(0, haloSigma)
        velocities[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo, 1] = randomGenerator.normal(0, haloSigma)
        velocities[nbDisk + nbBulge: nbDisk + nbBulge + nbHalo, 2] = randomGenerator.normal(0, haloSigma)
    return toParticleGroup(positions, velocities, masses, parameters.device)


def generate(parameters: SimulatorParameters) -> engine.ParticleGroup:
    if parameters.distributionType.upper() == "GALAXY":
        return generateGalaxyDistribution(parameters)
    return generateBasicDistribution(parameters)
