from sources.common.functions import *
from sources.galaxies import *


class ClusterEngine2D:
    def __init__(self, clusters=None, theta=1):
        self.gravitationCst = gravitationalConstant()
        if clusters is None:
            clusters = []
        self.clusters = clusters
        self.theta = theta
        # Massless Particles
        self.masslessSimulation = False
        self.masslessBool = (MasslessGalaxy, RingsMasslessGalaxy)
        self.masslessX, self.masslessV = None, None
        self.nbMassless = 0
        # Mass Clusters
        self.massesX, self.massesV, self.massesM = None, None, None
        self.nbMasses, self.nbCenters = 0, 0
        # Engine Variables
        self.quadTree = None
        # Load Objects
        self.loadObjects()

    def loadObjects(self):
        # Massive Clusters
        massesX = [cluster.positions for cluster in self.clusters if not isinstance(cluster, self.masslessBool)]
        massesV = [cluster.velocities for cluster in self.clusters if not isinstance(cluster, self.masslessBool)]
        massesM = [cluster.masses for cluster in self.clusters if not isinstance(cluster, self.masslessBool)]
        # Massless Clusters
        masslessObjects = [isinstance(cluster, self.masslessBool) for cluster in self.clusters]
        if True in masslessObjects:
            self.masslessSimulation = True
            # Loading Massless Particles
            masslessX = [cluster.positions for cluster in self.clusters if isinstance(cluster, self.masslessBool)]
            masslessV = [cluster.velocities for cluster in self.clusters if isinstance(cluster, self.masslessBool)]
            self.masslessX = np.concatenate(masslessX)
            self.masslessV = np.concatenate(masslessV)
            self.nbMassless = self.masslessX.shape[0]
            # Load Massless Cluster Centers
            masslessCentersX = [cluster.centerPosition for cluster in self.clusters]
            masslessCentersV = [cluster.centerVelocity for cluster in self.clusters]
            masslessCenterMasses = [cluster.centralMass for cluster in self.clusters]
            self.nbCenters = len(masslessCentersX)
            for i in range(self.nbCenters):
                massesX.append(masslessCentersX[i])
                massesV.append(masslessCentersV[i])
                massesM.append(masslessCenterMasses[i])
        self.massesX = np.array(massesX)
        self.massesV = np.array(massesV)
        self.massesM = np.array(massesM)
        self.nbMasses = self.massesX.shape[0]

    def acceleration(self, massesX, massesM, masslessX=None):
        quadTree = QuadTree(massesX, massesM, self.theta)
        # Masses particles
        massesA = np.zeros(masslessX.shape)
        n = massesX.shape[0]
        for j in range(n):
            massesA[j, :] = quadTree.forces(massesX[j, :])
        # Massless particles
        if self.masslessSimulation:
            masslessA = np.zeros(masslessX.shape)
            n = masslessX.shape[0]
            for j in range(n):
                masslessA[j, :] = quadTree.forces(masslessX[j, :])
            return massesA, masslessA
        else:
            return massesA

    def massesAcceleration(self, massesX, massesM):
        quadTree = QuadTree(massesX, massesM, self.theta)
        a = np.zeros(massesX.shape)
        n = massesX.shape[0]
        for j in range(n):
            a[j, :] = quadTree.forces(massesX[j, :])
        return a


class OctTree:
    def __init__(self, massesX, massesM, theta, limits=None):
        if limits is None:
            xMax = np.max(np.abs(massesX[:, 0])) * 1.01
            yMax = np.max(np.abs(massesX[:, 1])) * 1.01
            zMax = np.max(np.abs(massesX[:, 2])) * 1.01
            limits = (-xMax, xMax, -yMax, yMax, -zMax, zMax)
        if not type(massesX).__module__ == np.__name__:
            massesX = np.array([massesX])
        if not type(massesM).__module__ == np.__name__:
            massesM = np.array([massesM])
        self.massesX = massesX
        self.massesM = massesM
        self.gravitationCst = gravitationalConstant()
        self.nbObjects = self.massesX.shape[0]
        self.children = [None, None, None, None, None, None, None, None]
        self.x1, self.x2 = limits[0], limits[1]
        self.y1, self.y2 = limits[2], limits[3]
        self.z1, self.z2 = limits[4], limits[5]
        self.centreMassX = None
        self.centreMass = None
        self.theta = theta
        self.getCenter()
        self.buildChildren()

    def getCenter(self):
        xc = np.sum(self.massesX[:, 0])
        yc = np.sum(self.massesX[:, 1])
        zc = np.sum(self.massesX[:, 2])
        mc = np.sum(self.massesM)
        self.centreMass = mc
        self.centreMassX = np.array([xc, yc, zc])

    def buildChildren(self):
        # Top : towards +Z, North : towards +Y, East : towards +X
        # Each list :  contains positions [0] and masses [1]
        TNW, TNE, TSW, TSE = [[], []], [[], []], [[], []], [[], []]
        BNW, BNE, BSW, BSE = [[], [ ]], [[], []], [[], []], [[], []]
        # Calculating middle points
        mx = (self.x1 + self.x2) / 2
        my = (self.y1 + self.y2) / 2
        mz = (self.z1 + self.z2) / 2
        # Filling children recursively
        for i in range(self.nbObjects):
            if self.massesX[i, 0] > mx:  # East
                if self.massesX[i, 1] > my:  # North
                    if self.massesX[i, 2] > mz:  # Top
                        TNE[0].append(self.massesX[i, :])
                        TNE[1].append(self.massesM[i])
                    else:  # Bottom
                        BNE[0].append(self.massesX[i, :])
                        BNE[1].append(self.massesM[i])
                else:  # South
                    if self.massesX[i, 2] > mz:  # Top
                        TSE[0].append(self.massesX[i, :])
                        TSE[1].append(self.massesM[i])
                    else:  # Bottom
                        BSE[0].append(self.massesX[i, :])
                        BSE[1].append(self.massesM[i])
            else:  # West
                if self.massesX[i, 1] > my:  # North
                    if self.massesX[i, 2] > mz:  # Top
                        TNW[0].append(self.massesX[i, :])
                        TNW[1].append(self.massesM[i])
                    else:  # Bottom
                        BNW[0].append(self.massesX[i, :])
                        BNW[1].append(self.massesM[i])
                else:  # South
                    if self.massesX[i, 2] > mz:  # Top
                        TSW[0].append(self.massesX[i, :])
                        TSW[1].append(self.massesM[i])
                    else:  # Bottom
                        BSW[0].append(self.massesX[i, :])
                        BSW[1].append(self.massesM[i])
        if len(TNE[0]) > 0:
            self.children[0] = OctTree(TNE[0], TNE[1], self.theta, (mx, self.x2, my, self.y2, mz, self.z2))
        if len(BNE[0]) > 0:
            self.children[1] = OctTree(BNE[0], BNE[1], self.theta, (mx, self.x2, my, self.y2, self.z1, mz))
        if len(TSE[0]) > 0:
            self.children[2] = OctTree(TSE[0], TSE[1], self.theta, (mx, self.x2, self.y1, my, mz, self.z2))
        if len(BSE[0]) > 0:
            self.children[3] = OctTree(BSE[0], BSE[1], self.theta, (mx, self.x2, self.y1, my, self.z1, mz))
        if len(TNW[0]) > 0:
            self.children[4] = OctTree(TNW[0], TNW[1], self.theta, (self.x1, mx, my, self.y2, mz, self.z2))
        if len(BNW[0]) > 0:
            self.children[5] = OctTree(BNW[0], BNW[1], self.theta, (self.x1, mx, my, self.y2, self.z1, mz))
        if len(TSW[0]) > 0:
            self.children[6] = OctTree(TSW[0], TSW[1], self.theta, (self.x1, mx, self.y1, my, mz, self.z2))
        if len(BSW[0]) > 0:
            self.children[7] = OctTree(BSW[0], BSW[1], self.theta, (self.x1, mx, self.y1, my, self.z1, mz))

    def forces(self, position):
        if not type(position).__module__ == np.__name__:
            position = np.array([position])
        width = abs(self.x1 - self.x2)
        vector = self.centreMassX - position
        distance = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
        if distance != 0:
            ratio = width / distance
            if ratio < self.theta or self.allNoneChildren():
                acceleration = (self.gravitationCst * self.centreMass * vector) / (distance ** 3)
                return acceleration
            else:
                acceleration = np.array([0, 0, 0])
                for i in self.children:
                    if isinstance(i, OctTree):
                        acceleration += i.forces(position)
        else:
            return np.array([0, 0, 0])

    def allNoneChildren(self):
        for i in self.children:
            if i is not None:
                return False
        return True


class QuadTree:
    def __init__(self, massesX, massesM, theta, limits=None):
        if limits is None:
            xMax = np.max(np.abs(massesX[:, 0])) * 1.01
            yMax = np.max(np.abs(massesX[:, 1])) * 1.01
            limits = (-xMax, xMax, -yMax, yMax)
        if not type(massesX).__module__ == np.__name__:
            massesX = np.array([massesX])
        if not type(massesM).__module__ == np.__name__:
            massesM = np.array([massesM])
        self.massesX = massesX
        self.massesM = massesM
        self.gravitationCst = gravitationalConstant()
        self.nbObjects = self.massesX.shape[0]
        self.children = [None, None, None, None]
        self.x1, self.x2, self.y1, self.y2 = limits[0], limits[1], limits[2], limits[3]
        self.centreMassX = None
        self.centreMass = None
        self.theta = theta
        self.getCenter()
        self.buildChildren()

    def getCenter(self):
        xc = np.sum(self.massesX[:, 0])
        yc = np.sum(self.massesX[:, 1])
        mc = np.sum(self.massesM)
        self.centreMass = mc
        self.centreMassX = np.array([xc, yc])

    def buildChildren(self):
        # North : towards +Y, East : towards +X
        # Each list :  contains positions [0] and masses [1]
        NW, NE, SW, SE = [[], []], [[], []], [[], []], [[], []]
        # Calculating middle points
        mx = (self.x1 + self.x2) / 2
        my = (self.y1 + self.y2) / 2
        # Filling children recursively
        for i in range(self.nbObjects):
            if self.massesX[i, 0] > mx:  # East
                if self.massesX[i, 1] > my:  # North
                    NE[0].append(self.massesX[i, :])
                    NE[1].append(self.massesM[i])
                else:  # South
                    SE[0].append(self.massesX[i, :])
                    SE[1].append(self.massesM[i])
            else:  # West
                if self.massesX[i, 1] > my:  # North
                    NW[0].append(self.massesX[i, :])
                    NW[1].append(self.massesM[i])
                else:  # South
                    SW[0].append(self.massesX[i, :])
                    SW[1].append(self.massesM[i])
        if len(NE[0]) > 0:
            self.children[0] = QuadTree(NE[0], NE[1], self.theta, (mx, self.x2, my, self.y2))
        if len(SE[0]) > 0:
            self.children[1] = QuadTree(SE[0], SE[1], self.theta, (mx, self.x2, self.y1, my))
        if len(NW[0]) > 0:
            self.children[2] = QuadTree(NW[0], NW[1], self.theta, (self.x1, mx, my, self.y2))
        if len(SW[0]) > 0:
            self.children[3] = QuadTree(SW[0], SW[1], self.theta, (self.x1, mx, self.y1, my))

    def forces(self, position):
        if not type(position).__module__ == np.__name__:
            position = np.array([position])
        width = abs(self.x1 - self.x2)
        vector = position - self.centreMassX
        distance = np.sqrt(vector[0] ** 2 + vector[1] ** 2)
        if distance != 0:
            ratio = width / distance
            if ratio < self.theta or self.allNoneChildren():
                acceleration = (self.gravitationCst * self.centreMass * vector) / (distance ** 3)
                return acceleration
            else:
                acceleration = np.array([0, 0])
                for i in self.children:
                    if isinstance(i, QuadTree):
                        acceleration += i.forces(position)
        else:
            return np.array([0, 0, 0])

    def allNoneChildren(self):
        for i in self.children:
            if i is not None:
                return False
        return True
