import numpy as np

from sources.galaxies import *
import sys
import threading


class OctTree:
    def __init__(self, theta, limits=None, softLength=1/10**4):
        self.children = [None, None, None, None, None, None, None, None]
        self.massesX = None
        self.massesM = None
        self.gravitationCst = gravitationalConstant()
        self.nbObjects = 0
        self.x1, self.x2 = limits[0], limits[1]
        self.y1, self.y2 = limits[2], limits[3]
        self.z1, self.z2 = limits[4], limits[5]
        self.theta = theta
        self.softLength = softLength
        # Mass Center
        self.centreMassX = 0
        self.centreMass = np.array([0, 0, 0])

    def getCenter(self):
        Mc = np.sum(self.massesM)
        self.centreMass = Mc
        if Mc != 0:
            Rc = np.array([0.0, 0.0, 0.0])
            for i in range(self.nbObjects):
                Rc += self.massesM[i] * self.massesX[i, :]
            self.centreMassX = Rc / Mc

    def insert(self, position, mass):
        # Top : towards +Z, North : towards +Y, East : towards +X
        # Each list :  contains positions [0] and masses [1]
        # Calculating middle points
        mx = (self.x1 + self.x2) / 2
        my = (self.y1 + self.y2) / 2
        mz = (self.z1 + self.z2) / 2
        # Filling children recursively
        ###### ADDING PARTICLE TO SELF ######
        if self.nbObjects == 0:
            self.massesM = np.array([mass])
            self.massesX = np.array([position])
        else:
            self.children = [OctTree(self.theta, (mx, self.x2, my, self.y2, mz, self.z2)),  # TNE
                             OctTree(self.theta, (mx, self.x2, my, self.y2, self.z1, mz)),  # BNE
                             OctTree(self.theta, (mx, self.x2, self.y1, my, mz, self.z2)),  # TSE
                             OctTree(self.theta, (mx, self.x2, self.y1, my, self.z1, mz)),  # BSE
                             OctTree(self.theta, (self.x1, mx, my, self.y2, mz, self.z2)),  # BNW
                             OctTree(self.theta, (self.x1, mx, my, self.y2, self.z1, mz)),  # BNW
                             OctTree(self.theta, (self.x1, mx, self.y1, my, mz, self.z2)),  # BSW
                             OctTree(self.theta, (self.x1, mx, self.y1, my, self.z1, mz))]  # BSW
            self.massesM = np.append(self.massesM, np.array([mass]), 0)
            self.massesX = np.vstack([self.massesX, position])
            ###### ADDING PARTICLE TO CHILDREN ######
            if position[0] > mx:  # East
                if position[1] > my:  # North
                    if position[2] > mz:  # Top
                        self.children[0].insert(position, mass)
                    else:                 # Bottom
                        self.children[1].insert(position, mass)
                else:                 # South
                    if position[2] > mz:  # Top
                        self.children[2].insert(position, mass)
                    else:                 # Bottom
                        self.children[3].insert(position, mass)
            else:                 # West
                if position[1] > my:  # North
                    if position[2] > mz:  # Top
                        self.children[4].insert(position, mass)
                    else:  # Bottom
                        self.children[5].insert(position, mass)
                else:                 # South
                    if position[2] > mz:  # Top
                        self.children[6].insert(position, mass)
                    else:  # Bottom
                        self.children[7].insert(position, mass)
        ###### RECALCULATING CENTRE OF MASS ######
        self.nbObjects += 1
        self.getCenter()

    def forces(self, position):
        if not type(position).__module__ == np.__name__:
            position = np.array([position])
        width = abs(self.x1 - self.x2)
        vector = self.centreMassX - position
        distance = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
        if distance != 0:
            ratio = width / distance
            if ratio < self.theta or self.allNoneChildren():
                acceleration = (self.gravitationCst * self.centreMass * vector) / (distance + self.softLength) ** 3
                return acceleration
            else:
                acceleration = np.array([0.0, 0.0, 0.0])
                for i in self.children:
                    if isinstance(i, OctTree) and i.nbObjects != 0:
                        acceleration += i.forces(position)
        else:
            return np.array([0, 0, 0])

    def allNoneChildren(self):
        for i in self.children:
            if i is not None:
                return False
        return True


class QuadTree:
    def __init__(self, theta, limits=None, softLength=1/10**4):
        self.children = [None, None, None, None]
        self.massesX = None
        self.massesM = None
        self.gravitationCst = gravitationalConstant()
        self.nbObjects = 0
        self.x1, self.x2, self.y1, self.y2 = limits[0], limits[1], limits[2], limits[3]
        self.theta = theta
        self.softLength = softLength
        # Mass Center
        self.centreMassX = 0
        self.centreMass = np.array([0, 0])

    def getCenter(self):
        Mc = np.sum(self.massesM)
        self.centreMass = Mc
        if Mc != 0:
            Rc = np.array([0.0, 0.0])
            for i in range(self.nbObjects):
                Rc += self.massesM[i] * self.massesX[i, :]
            self.centreMassX = Rc / Mc

    def insert(self, position, mass):
        # North : towards +Y, East : towards +X
        # Each list :  contains positions [0] and masses [1]
        # Calculating middle points
        mx = (self.x1 + self.x2) / 2
        my = (self.y1 + self.y2) / 2
        ###### ADDING PARTICLE TO SELF ######
        if self.nbObjects == 0:
            self.massesM = np.array([mass])
            self.massesX = np.array([position])
        else:
            self.children = [QuadTree(self.theta, (mx, self.x2, my, self.y2), self.softLength),  # NE
                             QuadTree(self.theta, (mx, self.x2, self.y1, my), self.softLength),  # SE
                             QuadTree(self.theta, (self.x1, mx, my, self.y2), self.softLength),  # NW
                             QuadTree(self.theta, (self.x1, mx, self.y1, my), self.softLength)]  # SW
            self.massesM = np.append(self.massesM, np.array([mass]), 0)
            self.massesX = np.vstack([self.massesX, position])
            ###### ADDING PARTICLE TO CHILDREN ######
            if position[0] > mx:  # East
                if position[1] > my:  # North
                    self.children[0].insert(position, mass)
                else:                 # South
                    self.children[1].insert(position, mass)
            else:                 # West
                if position[1] > my:  # North
                    self.children[2].insert(position, mass)
                else:                 # South
                    self.children[3].insert(position, mass)
        ###### RECALCULATING CENTRE OF MASS ######
        self.nbObjects += 1
        self.getCenter()

    def forces(self, position):
        if not type(position).__module__ == np.__name__:
            position = np.array([position])
        width = abs(self.x1 - self.x2)
        vector = position - self.centreMassX
        distance = np.sqrt(vector[0] ** 2 + vector[1] ** 2)
        if distance != 0:
            ratio = width / distance
            if ratio < self.theta or self.allNoneChildren():
                acceleration = (self.gravitationCst * self.centreMass * vector) / (distance + self.softLength) ** 3
                return acceleration
            else:
                acceleration = np.array([0.0, 0.0])
                for i in self.children:
                    if isinstance(i, QuadTree) and i.nbObjects != 0:
                        acceleration += i.forces(position)
        else:
            return np.array([0, 0])

    def allNoneChildren(self):
        for i in self.children:
            if i is not None:
                return False
        return True


class ClusterEngine:
    def __init__(self, clusters=None, dimension=2, mode='BARNES_HUT', softeningLength=1/10000, theta=1, percentage=50):
        self.gravitationCst = gravitationalConstant()
        self.mode = mode
        self.dimension = dimension
        # Mass Clusters
        self.massesX, self.massesV, self.massesM = None, None, None
        self.nbMasses = 0
        # Engine Variables
        self.quadTree = None
        self.theta = None
        self.percentage = None
        self.softLength = softeningLength
        # Loading Clusters
        self.clusters = None
        if type(clusters) == list:
            self.feedClusters(clusters)
        if mode == 'BARNES_HUT':
            self.theta = theta
            self.acceleration = self.accelerationBarnesHut
        elif mode == 'ALTERED':
            self.percentage = percentage
            self.acceleration = self.accelerationAltered

    def feedClusters(self, clusters):
        self.clusters = clusters
        self.loadObjects()

    def setTheta(self, theta=1):
        self.theta = theta

    def loadObjects(self):
        # Massive Clusters
        massesX = [cluster.positions for cluster in self.clusters]
        massesV = [cluster.velocities for cluster in self.clusters]
        massesM = [cluster.masses for cluster in self.clusters]
        # Assigning Massive Particles Arrays
        self.massesX = np.concatenate(massesX)
        self.massesV = np.concatenate(massesV)
        self.massesM = np.concatenate(massesM)
        self.nbMasses = self.massesX.shape[0]
        sys.setrecursionlimit(self.nbMasses * 10)

    def accelerationBarnesHut(self, massesX, massesM):
        """Defective"""
        xMax = np.max(np.abs(massesX[:, 0])) * 1.01
        yMax = np.max(np.abs(massesX[:, 1])) * 1.01
        if self.dimension == 2:
            limits = (-xMax, xMax, -yMax, yMax)
            bhTree = QuadTree(self.theta, limits, self.softLength)
        elif self.dimension == 3:
            zMax = np.max(np.abs(massesX[:, 2])) * 1.01
            limits = (-xMax, xMax, -yMax, yMax, -zMax, zMax)
            bhTree = OctTree(self.theta, limits, self.softLength)
        n = massesX.shape[0]
        for i in range(n):
            bhTree.insert(massesX[i, :], massesM[i])
        # Masses particles
        massesA = np.zeros(massesX.shape)
        for j in range(n):
            massesA[j, :] = bhTree.forces(massesX[j, :])
        return massesA

    def accelerationAltered(self, massesX, massesM):
        massesA = np.zeros(massesX.shape)
        n = massesX.shape[0]
        limit = int(n * self.percentage / 100)
        for i in range(n):
            tags = np.concatenate([np.ones(shape=(limit,)), np.zeros(shape=(n - limit,))])
            np.random.shuffle(tags)
            for j in range(n):
                if j != i and tags[j] == 1:
                    vector = massesX[j, :] - massesX[i, :]
                    distance = np.linalg.norm(vector)
                    if distance != 0:
                        massesA[i, :] += (self.gravitationCst * massesM[j] * vector) / (distance + self.softLength) ** 3
        return massesA

    def compute(self, dt, method='EULER_EXPLICIT'):
        if method == 'EULER_EXPLICIT':
            self.computeEulerExplicit(dt)
        if method == 'EULER_SEMI_IMPLICIT':
            self.computeEulerSemiImplicit(dt)
        if method == 'RUNGE_KUTTA':
            self.computeRungeKutta(dt)

    def computeEulerExplicit(self, dt):
        # Calculating new positions
        newV = self.massesV + dt * self.acceleration(self.massesX, self.massesM)
        newX = self.massesX + dt * self.massesV
        # Updating positions
        self.massesX = newX
        self.massesV = newV

    def computeEulerSemiImplicit(self, dt):
        # Calculating new positions
        A = self.acceleration(self.massesX, self.massesM)
        newV = self.massesV + dt * A
        newX = self.massesX + dt * newV
        print(np.linalg.norm(A))
        print(np.linalg.norm(self.massesV))
        # Updating positions
        self.massesX = newX
        self.massesV = newV

    def computeRungeKutta(self, dt):
        # Calculating new positions
        k1_X = self.massesV * dt
        k1_V = self.acceleration(self.massesX, self.massesM) * dt
        k2_X = (self.massesV + k1_V / 2) * dt
        k2_V = self.acceleration(self.massesX + k1_X / 2, self.massesM) * dt
        k3_X = (self.massesV + k2_V / 2) * dt
        k3_V = self.acceleration(self.massesX + k2_X / 2, self.massesM) * dt
        k4_X = (self.massesV + k3_V) * dt
        k4_V = self.acceleration(self.massesX + k3_X, self.massesM) * dt
        # Updating positions
        newX = self.massesX + (k1_X + 2 * k2_X + 2 * k3_X + k4_X) / 6
        newV = self.massesV + (k1_V + 2 * k2_V + 2 * k3_V + k4_V) / 6
        self.massesX = newX
        self.massesV = newV