import os
import imageio
import time

import numpy as np
import matplotlib.pyplot as plt

from sources.common.constants import *


class ObjectCluster:
    def __init__(self, positions, velocities, darkPercentage):
        if not type(positions).__module__ == np.__name__:
            positions = np.array([positions])
        if not type(velocities).__module__ == np.__name__:
            velocities = np.array([velocities])
        assert positions.shape == velocities.shape
        assert darkPercentage >= 0, 'Dark Matter Percentage needs to be between 0 and 100'
        assert darkPercentage <= 100, 'Dark Matter Percentage needs to be between 0 and 100'
        self.gravitationCst = gravitationalConstant()
        self.darkPercentage = darkPercentage
        self.positions = positions
        self.velocities = velocities
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
        self.positions += xc
        self.velocities += vc
        self.initialized = True


class DiskGalaxy2D(ObjectCluster):
    def __init__(self, nbStars, radius, darkPercentage):




        positions, velocities, masses = None, None, None
        super().__init__(positions, velocities)

