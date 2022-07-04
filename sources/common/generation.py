import numpy as np
import matplotlib.pyplot as plt

from sources.common.gravitation import *


def particleRing(nb, radius, gravityCst, mass):
    particles = []
    velocities = []
    theta = 0
    arcLength = (2 * np.pi) / nb
    v = np.sqrt(gravityCst * mass / radius)
    while len(particles) < nb:
        angle = theta * arcLength
        beta = angle + np.pi / 2
        theta += 1
        particles.append([radius * np.cos(angle), radius * np.sin(angle)])
        velocities.append([v * np.cos(beta), v * np.sin(beta)])
    return np.array(particles), np.array(velocities)


def generateDisk2D(nbStars, radius, mass, gravityCst, seed=None):
    np.random.seed(seed)

    # Calculating positions
    positions = np.zeros(shape=(nbStars, 2))
    distances = np.random.random((nbStars,)) * radius
    angles = np.random.random((nbStars,)) * 2 * np.pi
    positions[:, 0] = np.cos(angles) * distances
    positions[:, 1] = np.sin(angles) * distances

    # Calculating speeds
    velocities = np.zeros(shape=(nbStars, 2))
    masses = np.random.random((nbStars,)) * mass
    for i in range(nbStars):
        mask = distances < distances[i]
        internalMass = np.sum(masses[mask])
        escVelocityNorm = np.sqrt(gravityCst * internalMass / distances[i])
        velocities[i, 0] = escVelocityNorm * np.cos(angles[i] + np.pi / 2)
        velocities[i, 1] = escVelocityNorm * np.sin(angles[i] + np.pi / 2)
    return positions, velocities, masses


def generateDisk3D(nbStars, radius, mass, zOffsetMax, gravityCst, seed=None):
    np.random.seed(seed)

    # Calculating positions
    positions = np.zeros(shape=(nbStars, 3))
    distances = np.random.random((nbStars,))
    zOffsets = (np.random.random((nbStars,)) - 0.5) * 2 * zOffsetMax * (np.ones_like(distances) - np.sqrt(distances))
    distances = distances * radius
    angles = np.random.random((nbStars,)) * 2 * np.pi
    positions[:, 0] = np.cos(angles) * distances
    positions[:, 1] = np.sin(angles) * distances
    positions[:, 2] = zOffsets

    # Calculating speeds
    velocities = np.zeros(shape=(nbStars, 3))
    masses = np.random.random((nbStars,)) * mass
    for i in range(nbStars):
        mask = distances > distances[i]
        internalMass = np.sum(masses[mask])
        velNorm = np.sqrt(gravityCst * internalMass / distances[i])
        velocities[i, 0] = velNorm * np.cos(angles[i])
        velocities[i, 1] = velNorm * np.sin(angles[i])
        velocities[i, 2] = np.zeros_like(velocities[i, 2])
    return positions, velocities, masses


def generateArms2D(nbStars, nbArms, radius, armOffset, mass, rotFactor, gravityCst, seed=None):
    np.random.seed(seed)
    armSeparationDistance = 2 * np.pi / nbArms
    distances = np.random.random((nbStars,)) ** 2
    angles = np.random.random((nbStars,)) * 2 * np.pi

    # Calculating arm offsets
    armOffsets = np.random.random((nbStars,)) * armOffset
    armOffsets = armOffsets - armOffset / 2
    armOffsets = armOffsets / distances
    squaredArmOffsets = armOffsets ** 2
    mask = armOffsets < 0
    squaredArmOffsets[mask] = -1 * squaredArmOffsets[mask]
    armOffsets = squaredArmOffsets

    # Rotation angles
    rotations = distances * rotFactor
    for i in range(nbStars):
        angles[i] = int(angles[i] / armSeparationDistance)
    angles = angles * armSeparationDistance + armOffsets + rotations

    # Calculating positions
    positions = np.zeros(shape=(nbStars, 2))
    positions[:, 0] = np.cos(angles) * distances * radius
    positions[:, 1] = np.sin(angles) * distances * radius

    # Calculating speeds
    velocities = np.zeros(shape=(nbStars, 2))
    masses = np.random.random((nbStars,)) * mass
    for i in range(nbStars):
        mask = distances > distances[i]
        internalMass = np.sum(masses[mask])
        velNorm = np.sqrt(gravityCst * internalMass / distances[i])
        velocities[i, 0] = velNorm * np.cos(angles[i])
        velocities[i, 0] = velNorm * np.sin(angles[i])
    return positions, velocities, masses
