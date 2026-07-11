import numpy as np
from PyQt5.QtCore import QThread, pyqtSignal, QObject
from concurrent.futures import ThreadPoolExecutor, as_completed

from src.core.engine.integrators import EulerExplicit, RK4
from src.core.engine.calculators import BarnesHutCalculator, NewtonCalculator
from src.core.engine.initializers import Initializer
from src.core.engine.parameters import SimulatorParameters


class NBodySimulator(QObject):
    positionsReady = pyqtSignal(dict)
    simulationFinished = pyqtSignal()

    def __init__(self, parameters: SimulatorParameters):
        super().__init__()
        self.parameters = parameters
        self.timeStep = parameters.timeStep
        self.gravitationalConstant = parameters.gravitationalConstant
        self.is3D = self._getIs3D()
        self.initializer = Initializer(gravitationalConstant=self.gravitationalConstant)
        self.calculator = self._createCalculator(parameters)
        self.integrator = self._createIntegrator(parameters)
        self.positions = self.velocities = self.masses = None
        self.time = 0.0
        self.isRunning = False
        self._stopRequested = False

    def setParameters(self, parameters: SimulatorParameters):
        self.parameters = parameters
        self.gravitationalConstant = parameters.gravitationalConstant
        self.timeStep = parameters.timeStep
        self.is3D = self._getIs3D()
        self.initializer = Initializer(gravitationalConstant=self.gravitationalConstant)
        self.calculator = self._createCalculator(parameters)
        self.integrator = self._createIntegrator(parameters)
        self.positions = self.velocities = self.masses = None

    def _getIs3D(self):
        if self.parameters.distributionType.lower() == "basic":
            return self.parameters.basicDistributionParameters.is3D
        else:
            return self.parameters.galaxyDistributionParameters.is3D

    @staticmethod
    def _createIntegrator(parameters: SimulatorParameters):
        integratorType = parameters.integratorType.upper()
        if integratorType == "EULER_EXPLICIT":
            return EulerExplicit(parameters.timeStep)
        else:
            return RK4(parameters.timeStep)

    @staticmethod
    def _createCalculator(parameters: SimulatorParameters):
        calculatorType = parameters.calculatorType.upper()
        if calculatorType == "BARNES_HUT":
            return BarnesHutCalculator(theta=parameters.theta, gravitationalConstant=parameters.gravitationalConstant, is3D=parameters.is3D)
        else:
            return NewtonCalculator(gravitationalConstant=parameters.gravitationalConstant, is3D=parameters.is3D)

    def initialize(self):
        distributionType = self.parameters.distributionType.lower()
        if distributionType == "basic":
            p = self.parameters.basicDistributionParameters
            self.positions, self.velocities, self.masses = self.initializer.basicDistribution(numParticles=p.numParticles, positionScale=p.positionScale, velocityScale=p.velocityScale, is3D=p.is3D)
        elif distributionType == "galaxy":
            p = self.parameters.galaxyDistributionParameters
            self.positions, self.velocities, self.masses = self.initializer.galaxyDistribution(numParticles=p.numParticles, totalMass=p.totalMass, radius=p.radius, height=p.height, is3D=p.is3D)
        else:
            raise ValueError("Unsupported distributionType")
        print(f"Initialized {distributionType} simulation with {self.parameters.numParticles if hasattr(self.parameters, 'numParticles') else 'N'} particles.")

    def run(self):
        print('NBodySimulator thread started')
        self.isRunning = True
        self._stopRequested = False
        try:
            self.initialize()
            stepCount = 0
            while self.isRunning and not self._stopRequested:
                print(stepCount)
                self.positions, self.velocities = self.integrator.step(self.positions, self.velocities, self.masses, self.calculator)
                self.time += self.timeStep
                groups = {"main": self.positions.copy()}
                self.positionsReady.emit(groups)
                stepCount += 1
        except Exception as e:
            print("Simulation error:", e)
        finally:
            self.isRunning = False
            self.simulationFinished.emit()

    def stop(self):
        self._stopRequested = True
        self.isRunning = False

    def step(self):
        if self.positions is None:
            self.initialize()
        self.positions, self.velocities = self.integrator.step(
            self.positions, self.velocities, self.masses, self.calculator)
        self.time += self.timeStep
        return self.positions, self.velocities
