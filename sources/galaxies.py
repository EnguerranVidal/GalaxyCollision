import os
import imageio
import time

import numpy as np
import matplotlib.pyplot as plt

from sources.common.functions import *


class ObjectCluster:
    def __init__(self, positions, velocities):
        if not type(positions).__module__ == np.__name__:
            positions = np.array([positions])
        if not type(velocities).__module__ == np.__name__:
            velocities = np.array([velocities])
        assert positions.shape == velocities.shape
        self.gravitationCst = gravitationalConstant()
        self.positions = positions
        self.velocities = velocities
        self.nbParticles = self.positions.shape[0]
        self.initialized = False

    def initialState(self, xc, vc):
        if not type(xc).__module__ == np.__name__:
            xc = np.array([xc])
        if not type(vc).__module__ == np.__name__:
            vc = np.array([vc])
        self.positions += xc
        self.velocities += vc
        self.initialized = True


class MasslessGalaxy(ObjectCluster):
    def __init__(self, positions, velocities, centralMass, haloRadius):
        super(MasslessGalaxy, self).__init__(positions, velocities)
        self.center_position = None
        self.center_velocity = None
        self.centralMass = centralMass
        self.haloRadius = haloRadius
        self.haloVelocity = 2 * np.sqrt(self.gravitationCst * self.centralMass / self.haloRadius)

    def initialState(self, xc, vc):
        if not type(xc).__module__ == np.__name__:
            xc = np.array([xc])
        if not type(vc).__module__ == np.__name__:
            vc = np.array([vc])
        self.positions += xc
        self.velocities += vc
        self.center_position = xc
        self.center_velocity = vc
        self.initialized = True

    def interiorMass(self, radii):
        indices = radii < self.haloRadius
        masses = np.zeros_like(radii)
        if masses[indices].shape != (0,):
            masses[indices] = self.haloVelocity ** 2 * radii[indices] ** 3 / \
                              (self.gravitationCst * (radii[indices] + self.haloRadius) ** 2)
        if masses[~indices].shape != (0,):
            masses[~indices] = self.centralMass
        return masses

    def density(self, radii):
        radii_in = 0.99 * radii
        radii_out = 1.01 * radii
        M_in = self.interiorMass(radii_in)
        M_out = self.interiorMass(radii_out)
        dM = M_out - M_in
        dV = (4 / 3) * np.pi * (radii_out ** 3 - radii_in ** 3)
        return dM / dV

    def dynamicFriction(self, radii, velocities, centralVelocity, mass):
        rho = self.density(radii)
        velocityVec = velocities - centralVelocity
        velocities = np.linalg.norm(velocityVec)
        return -4 * np.pi * self.gravitationCst ** 2 * 3 * mass * rho * velocityVec / (1 + velocities) ** 3


class RingsMasslessGalaxy(MasslessGalaxy):
    def __init__(self, radii, particles, centralMass, haloRadius):
        assert len(radii) == len(particles), "len(radii): " + str(len(radii)) + " not equal to len(particles): " + str(
            len(particles))
        G = gravitationalConstant()
        positions = []
        velocities = []
        for i in range(len(particles)):
            X, V = particleRing(particles[i], radii[i], G, centralMass)
            positions.append(X)
            velocities.append(V)
        positions = np.concatenate(positions)
        velocities = np.concatenate(velocities)
        super(RingsMasslessGalaxy, self).__init__(positions, velocities, centralMass, haloRadius)