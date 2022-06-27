import os
import imageio
import time

import numpy as np
import matplotlib.pyplot as plt

from sources.common.gravitation import gravitationalConstant
from sources.common.generation import generateDisk2D


class ObjectCluster2D:
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
        self.initPosition = None
        self.initVelocity = None
        self.initialState([0, 0], [0, 0])
        self.masses = masses
        self.nbParticles = self.positions.shape[0]
        self.initialized = False

    def initialState(self, xc, vc):
        if isinstance(xc, list):
            xc = np.array([xc])
            xc = xc.reshape(2, )
        if isinstance(vc, list):
            vc = np.array([vc])
            vc = vc.reshape(2, )
        xc = np.array([xc])
        vc = np.array([vc])
        self.initPosition = xc
        self.initVelocity = vc
        self.positions += xc
        self.velocities += vc
        self.initialized = True


class DiskGalaxy2D(ObjectCluster2D):
    def __init__(self, name, nbStars, radius, darkPercentage, mass=10, seed=None):
        gravConst = gravitationalConstant()
        positions, velocities, masses = generateDisk2D(nbStars, radius, mass, gravConst, seed)
        super().__init__(name, positions, velocities, masses, darkPercentage)

