from sources.galaxies import *
from sources.engines import *
from sources.simulators import *
from sources.common.graphs import *


try:
    import matplotlib
    matplotlib.use("TkAgg")
    import cupy as np
    print('CUPY VERSION DETECTED')
except ImportError:
    import numpy as np


def testGeneration():
    galaxy = DiskGalaxy3D('MILKY_WAY', 10000)
    plotCumulativeMass(galaxy)


def testUniform():
    plotUniformSphere()


def singleGalaxyTest():
    execDir = os.path.dirname(os.path.abspath(__file__))
    clusters = [DiskGalaxy2D('MILKY_WAY', 1000, 1, 75, mass=1)]
    sim = Simulator(clusters, execDir, mode='BARNES_HUT')
    sim.run(0.1, 50, method='EULER_SEMI_IMPLICIT')
    createGIF2D(sim.sessionPath, size=5)


def testCreationGIF():
    path = 'D:\\_PROJECTS\\PYTHON PROJECTS\\Galaxy_Merger\\logs\\20220630224442'
    createGIF2D(path, size=(100, 90))


if __name__ == "__main__":
    testGeneration()
