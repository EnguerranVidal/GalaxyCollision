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
        self.nbClusters = len(self.clusters)
        self.tags = np.concatenate([np.full((self.clusters[i].nbParticles,), i) for i in range(len(self.clusters))])
        self.simulatorDimension = max([len(i.initPosition) for i in self.clusters])
        self.simulatorTime = 0

        ###### Data Save System ######
        self.isNew = True
        self.nbSaves = 0
        self.sessionPath = None
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
            ###### POSITIONS FILE ######
            self.snapshotsFile = os.path.join(self.sessionPath, 'ClustersPositions.txt')
            self.saveState()
            ###### CLUSTERS FILE ######
            clustersLines = ['CLUSTER\tNB_PARTICLES\tDARK_MATTER\tINITIAL_X\tINITIAL_V\tDIMENSION\tNB_FRAMES\n']
            for i in range(self.nbClusters):
                name = self.clusters[i].name
                nbParticles = str(self.clusters[i].nbParticles)
                darkPercentage = str(self.clusters[i].darkPercentage)
                initialX, initialV = '', ''
                dimension = len(self.clusters[i].initPosition)
                for j in range(self.simulatorDimension):
                    initialX += str(self.clusters[i].initPosition[j]) + ' '
                    initialV += str(self.clusters[i].initVelocity[j]) + ' '
                # Creating Line
                line = name + '\t' + nbParticles + '\t' + darkPercentage + '\t' + initialX[:-1]
                line += '\t' + initialV[:-1] + '\t' + str(dimension) + str(self.nbSaves) + '\n'
                clustersLines.append(line)
            self.clustersFile = os.path.join(self.sessionPath, 'ClustersConfiguration.txt')
            with open(self.clustersFile, 'w') as file:
                for i in range(len(clustersLines)):
                    file.write(clustersLines[i])

    def saveState(self):
        ###### SAVE PARTICLES' POSITIONS ######
        modes = ['a', 'w']
        X = self.engine.massesX[:, 0].tolist()
        Y = self.engine.massesX[:, 1].tolist()
        Vx = self.engine.massesV[:, 0].tolist()
        Vy = self.engine.massesV[:, 1].tolist()
        with open(self.snapshotsFile, modes[int(self.isNew)]) as file:
            file.write(str(self.simulatorTime) + '\n')
            file.write(' '.join(X) + '\n')
            file.write(' '.join(Y) + '\n')
            if self.simulatorDimension == 3:
                Z = self.engine.massesX[:, 2].tolist()
                file.write(' '.join(Z) + '\n')
            file.write(' '.join(Vx) + '\n')
            file.write(' '.join(Vy) + '\n')
            if self.simulatorDimension == 3:
                Vz = self.engine.massesV[:, 2].tolist()
                file.write(' '.join(Vz) + '\n')
        ###### CHANGE NUMBER OF SAVES ######
        self.nbSaves += 1
        with open(self.clustersFile, 'r') as file:
            lines = file.readlines()
        for i in range(1, self.nbClusters + 1):
            line = lines[i].strip('\n')
            line = line.split('\t')
            line[-1] = str(int(line[-1]) + 1)
            line = '\t'.join(line) + '\n'
            lines[i] = line
        with open(self.clustersFile, 'w') as file:
            for i in range(len(lines)):
                file.write(lines[i])

    def run(self, dt=0.1, T=10, method='EULER_SEMI_IMPLICIT'):
        if self.isNew:
            self.initializeLogs()
        initial_time = self.simulatorTime
        print("########## Beginning Calculations ##########")
        while self.simulatorTime < initial_time + T:
            self.engine.compute(dt, method=method)
            self.simulatorTime = self.simulatorTime + dt
            self.saveState()
