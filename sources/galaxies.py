import os
import imageio
import time

import numpy as np
import matplotlib.pyplot as plt

from sources.common.gravitation import gravitationalConstant
from sources.common.generation import generateDisk2D, generateUniformSphere, generateDisk3D


class ObjectCluster:
    def __init__(self, name, positions, velocities, masses, darkPercentage):
        if not type(positions).__module__ == np.__name__:
            positions = np.array([positions])
        if not type(velocities).__module__ == np.__name__:
            velocities = np.array([velocities])
        if not type(masses).__module__ == np.__name__:
            masses = np.array([masses])
        assert positions.shape == velocities.shape
        assert darkPercentage >= 0, 'Dark Matter Percentage needs to be between 0 and 100'
        assert darkPercentage <= 100, 'Dark Matter Percentage needs to be between 0 and 100'

        self.gravitationCst = gravitationalConstant()
        self.darkPercentage = darkPercentage
        self.name = name
        self.positions = positions
        self.velocities = velocities
        self.masses = masses

        self.dimension = self.positions.shape[1]
        self.nbParticles = self.positions.shape[0]
        self.initPosition = None
        self.initVelocity = None
        self.initialState([0] * self.dimension, [0] * self.dimension)
        self.initialized = False

    def initialState(self, xc, vc):
        if isinstance(xc, list):
            xc = np.array([xc])
            xc = xc.reshape(xc.shape[1], )
        if isinstance(vc, list):
            vc = np.array([vc])
            vc = vc.reshape(vc.shape[1], )
        xc = np.array([xc])
        vc = np.array([vc])
        self.initPosition = xc
        self.initVelocity = vc
        self.positions += xc
        self.velocities += vc
        self.initialized = True


class DiskGalaxy2D(ObjectCluster):
    def __init__(self, name, nbStars=1000, radius=1, darkPercentage=75, mass=1, seed=None):
        gravConst = gravitationalConstant()
        positions, velocities, masses = generateDisk2D(nbStars, radius, mass, gravConst, seed)
        super().__init__(name, positions, velocities, masses, darkPercentage)


class UniformSphericalCluster3D(ObjectCluster):
    def __init__(self, name, nbStars=1000, radius=1, darkPercentage=75, mass=1, seed=None):
        gravConst = gravitationalConstant()
        positions, velocities, masses = generateUniformSphere(nbStars, radius, mass, gravConst, seed)
        super().__init__(name, positions, velocities, masses, darkPercentage)


class DiskGalaxy3D(ObjectCluster):
    def __init__(self, name, nbStars=1000, radius=1, darkPercentage=75, mass=1, zOffsetMax=0.05, seed=None):
        gravConst = gravitationalConstant()
        positions, velocities, masses = generateDisk3D(nbStars, radius, mass, zOffsetMax, gravConst, seed)
        super().__init__(name, positions, velocities, masses, darkPercentage)

