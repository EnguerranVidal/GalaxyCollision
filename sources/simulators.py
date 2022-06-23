import os

import numpy as np
import matplotlib.pyplot as plt

from sources.common.other import sessionName
from sources.common.gravitation import gravitationalConstant
from sources.common.generation import generateDisk2D
from sources.galaxies import DiskGalaxy2D
from sources.engines import ClusterEngine2D, QuadTree


class Simulator:
    def __init__(self, clusters, engine, path):
        self.execDir = path
        self.engine = engine
        self.clusters = clusters
        self.tags = np.concatenate([np.full((self.clusters[i].nbParticles,), i) for i in range(len(self.clusters))])
        self.nbClusters = len(self.clusters)
        self.time = 0

        # Data Save System
        self.isNew = True
        self.sessionPath = None
        self.nSaves = 0

        # Files
        self.massesFile = None
        self.snapshotsFile = None
        self.clustersFile = None
        self.simulationFile = None

    def reset(self):
        self.__init__(self.clusters, self.engine, self.execDir)

    def initializeLogs(self):
        if self.isNew:
            logsPath = os.path.join(self.execDir + 'logs')
            if not os.path.exists(logsPath):
                os.mkdir(logsPath)
            logFile = sessionName()
            self.sessionPath = os.path.join(logsPath, logFile)
            if not os.path.exists(self.sessionPath):
                os.mkdir(self.sessionPath)

            ###### MASSES FILE ######
            massesLines = ['CLUSTER' + '\t' + 'MASSES' + '\n']
            for i in range(self.nbClusters):
                massesLines.append([])
                name = self.clusters[i].name
                for j in range(self.clusters[i].nbParticles):
                    massesLines[i].append(str(self.clusters[i].massesM[j]))
                massesLines[i] = name + '\t' + ' '.join(massesLines[i]) + '\n'
            self.massesFile = os.path.join(self.sessionPath, 'MassesConfiguration.txt')
            with open(self.massesFile, 'w') as file:
                for i in range(len(massesLines)):
                    file.write(massesLines[i])

            ###### CLUSTERS FILE ######
            clustersLines = ['CLUSTER' + '\t' + 'NB_PARTICLES' + '\t' + 'DARK_MATTER' +
                             '\t' 'INITIAL_X' + '\t' + 'INITIAL_V' + '\n']
            for i in range(self.nbClusters):
                name = self.clusters[i].name
                nbParticles = str(self.clusters[i].nbParticles)
                darkPercentage = str(self.clusters[i].darkPercentage)
                initialX, initialV = '', ''
                for j in range(len(self.clusters[i].initPosition)):
                    initialX += str(self.clusters[i].initPosition[j]) + ' '
                    initialV += str(self.clusters[i].initVelocity[j]) + ' '
                line = name + '\t' + nbParticles + '\t' + darkPercentage + '\t'
                line += initialX[:-1] + '\t' + initialV[:-1] + '\n'
                clustersLines.append(line)
            self.clustersFile = os.path.join(self.sessionPath, 'ClustersConfiguration.txt')
            with open(self.clustersFile, 'w') as file:
                for i in range(len(clustersLines)):
                    file.write(clustersLines[i])

    def saveState(self):
        self.nSaves += 1
        pass

    def run(self, dt=0.1, T=10, method='EULER_SEMI_IMPLICIT'):
        initial_time = self.time
        while self.time < initial_time + T:
            self.engine.compute(dt, method=method)
            self.time = self.time + dt

