from sources.galaxies import *
import sys
import threading


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
        BNW, BNE, BSW, BSE = [[], []], [[], []], [[], []], [[], []]
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
    def __init__(self, theta, limits=None):
        self.children = [None, None, None, None]
        self.massesX = None
        self.massesM = None
        self.gravitationCst = gravitationalConstant()
        self.nbObjects = 0
        self.x1, self.x2, self.y1, self.y2 = limits[0], limits[1], limits[2], limits[3]
        self.theta = theta
        # Mass Center
        self.centreMassX = None
        self.centreMass = None

    def getCenter(self):
        mc = np.sum(self.massesM)
        self.centreMass = mc
        xc = np.sum(self.massesX[:, 0])
        yc = np.sum(self.massesX[:, 1])
        self.centreMassX = np.array([xc, yc])

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
            self.children = [QuadTree(self.theta, (mx, self.x2, my, self.y2)),  # NE
                             QuadTree(self.theta, (mx, self.x2, self.y1, my)),  # SE
                             QuadTree(self.theta, (self.x1, mx, my, self.y2)),  # NW
                             QuadTree(self.theta, (self.x1, mx, self.y1, my))]  # SW
            np.append(self.massesM, mass, 0)
            np.append(self.massesX, position, 0)
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


class ClusterEngine2D:
    def __init__(self, clusters=None, mode='BARNES_HUT', theta=1, percentage=50):
        self.gravitationCst = gravitationalConstant()
        self.mode = mode
        # Mass Clusters
        self.massesX, self.massesV, self.massesM = None, None, None
        self.nbMasses = 0
        # Engine Variables
        self.quadTree = None
        self.theta = None
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
        limits = (-xMax, xMax, -yMax, yMax)
        quadTree = QuadTree(self.theta, limits)
        n = massesX.shape[0]
        for i in range(n):
            quadTree.insert(massesX[i, :], massesM[i])
        # Masses particles
        massesA = np.zeros(massesX.shape)
        for j in range(n):
            massesA[j, :] = quadTree.forces(massesX[j, :])
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
                        massesA[i, :] += (self.gravitationCst * massesM[j] * vector) / (distance ** 3)
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
        newV = self.massesV + dt * self.acceleration(self.massesX, self.massesM)
        newX = self.massesX + dt * newV
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
