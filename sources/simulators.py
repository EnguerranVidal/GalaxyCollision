import os
import time

import numpy as np
import matplotlib.pyplot as plt

from sources.common.other import sessionName, progressBar
from sources.common.gravitation import gravitationalConstant
from sources.common.generation import generateDisk2D
from sources.galaxies import DiskGalaxy2D
from sources.engines import ClusterEngine, QuadTree


class Simulator:
    def __init__(self, clusters, path, mode='BARNES_HUT',
                 softeningLength=1/10000, theta=1, percentage=50):
        self.execDir = path
        self.clusters = clusters
        self.nbClusters = len(self.clusters)
        self.tags = np.concatenate([np.full((self.clusters[i].nbParticles,), i) for i in range(len(self.clusters))])
        self.simulatorDimension = max([i.initPosition.shape[1] for i in self.clusters])
        self.engine = ClusterEngine(self.clusters, self.simulatorDimension, mode, softeningLength, theta, percentage)
        self.simulatorTime = 0

        ###### Data Save System ######
        self.isNew = True
        self.nbSaves = 0
        self.sessionPath = None
        self.massesFile = None
        self.snapshotsFile = None
        self.clustersFile = None
        self.simulationFile = None

    def reset(self, mode='BARNES_HUT', softeningLength=1/10000, theta=1, percentage=50):
        self.__init__(self.clusters, self.execDir, mode, softeningLength, theta, percentage)

    def initializeLogs(self):
        if self.isNew:
            logsPath = os.path.join(self.execDir, 'logs')
            if not os.path.exists(logsPath):
                os.mkdir(logsPath)
            logFile = sessionName()
            self.sessionPath = os.path.join(logsPath, logFile)
            if not os.path.exists(self.sessionPath):
                os.mkdir(self.sessionPath)
            ###### MASSES FILE ######
            massesLines = []
            for i in range(self.nbClusters):
                massesLines.append([])
                for j in range(self.clusters[i].nbParticles):
                    massesLines[i].append(str(self.clusters[i].masses[j]))
                massesLines[i] = ' '.join(massesLines[i]) + '\n'
            self.massesFile = os.path.join(self.sessionPath, 'MassesConfiguration.config')
            with open(self.massesFile, 'w') as file:
                for i in range(len(massesLines)):
                    file.write(massesLines[i])
            ###### CLUSTERS FILE ######
            clustersLines = ['CLUSTER NB_PARTICLES DARK_MATTER INITIAL_X INITIAL_V DIMENSION NB_FRAMES\n']
            for i in range(self.nbClusters):
                name = self.clusters[i].name
                nbParticles = str(self.clusters[i].nbParticles)
                darkPercentage = str(self.clusters[i].darkPercentage)
                initialX, initialV = '', ''
                dimension = len(self.clusters[i].initPosition[0])
                for j in range(dimension):
                    initialX += str(self.clusters[i].initPosition[0][j]) + ','
                    initialV += str(self.clusters[i].initVelocity[0][j]) + ','
                # Creating Line
                line = name + ' ' + nbParticles + ' ' + darkPercentage + ' ' + initialX[:-1]
                line += ' ' + initialV[:-1] + ' ' + str(dimension) + ' ' + str(self.nbSaves) + '\n'
                clustersLines.append(line)
            self.clustersFile = os.path.join(self.sessionPath, 'ClustersConfiguration.config')
            with open(self.clustersFile, 'w') as file:
                for i in range(len(clustersLines)):
                    file.write(clustersLines[i])
            ###### POSITIONS FILE ######
            self.snapshotsFile = os.path.join(self.sessionPath, 'ClustersPositions.config')
            self.saveState()

    def saveState(self):
        modes = ['a', 'w']
        ###### SAVE PARTICLES' POSITIONS ######
        X = self.engine.massesX[:, 0].tolist()
        Y = self.engine.massesX[:, 1].tolist()
        Vx = self.engine.massesV[:, 0].tolist()
        Vy = self.engine.massesV[:, 1].tolist()
        with open(self.snapshotsFile, modes[int(self.isNew)]) as file:
            file.write(str(self.simulatorTime) + '\n')
            file.write(' '.join(map(str, X)) + '\n')
            file.write(' '.join(map(str, Y)) + '\n')
            if self.simulatorDimension == 3:
                Z = self.engine.massesX[:, 2].tolist()
                file.write(' '.join(map(str, Z)) + '\n')
            file.write(' '.join(map(str, Vx)) + '\n')
            file.write(' '.join(map(str, Vy)) + '\n')
            if self.simulatorDimension == 3:
                Vz = self.engine.massesV[:, 2].tolist()
                file.write(' '.join(map(str, Vz)) + '\n')
        ###### CHANGE NUMBER OF SAVES ######
        self.nbSaves += 1
        with open(self.clustersFile, 'r') as file:
            lines = file.readlines()
        for i in range(1, self.nbClusters + 1):
            line = lines[i].strip('\n')
            line = line.split()
            line[-1] = str(int(line[-1]) + 1)
            line = ' '.join(line) + '\n'
            lines[i] = line
        with open(self.clustersFile, 'w') as file:
            for i in range(len(lines)):
                file.write(lines[i])

    def run(self, dt=0.1, T=10, method='EULER_SEMI_IMPLICIT'):
        if self.isNew:
            self.initializeLogs()
            self.isNew = not self.isNew
        print("########## Beginning Calculations ##########")
        # Progress Bar Displaying
        timeStamp0 = time.time()
        step, limit = 0, int(T / dt)
        progressBar(step, limit, 0, width=30)
        # MAIN LOOP
        for step in range(limit):
            # Computation
            self.engine.compute(dt, method=method)
            self.simulatorTime = self.simulatorTime + dt
            # Data Save
            self.saveState()
            # Progress Bar Update
            timeDelta = time.time() - timeStamp0
            progressBar(step + 1, limit, timeDelta, width=30)
        print("\n########## Calculations Finished ##########")
