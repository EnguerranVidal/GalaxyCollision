from PyQt5.QtCore import QThread, pyqtSignal, QObject, pyqtSlot

from src.core.engine.integrators import EulerExplicit, RK4
from src.core.engine.calculators import NewtonCalculator
from src.core.engine.distributions import Distribution, GalaxyDistribution
from src.core.engine.parameters import SimulatorParameters


class NBodySimulator(QObject):
    positionsReady = pyqtSignal(dict)
    simulationFinished = pyqtSignal()

    def __init__(self, parameters: SimulatorParameters):
        super().__init__()
        self.parameters = parameters
        self.timeStep = parameters.timeStep
        self.gravitationalConstant = parameters.gravitationalConstant
        self.calculator = self._createCalculator(parameters)
        self.integrator = self._createIntegrator(parameters)
        self.particles = self._createDistribution(parameters)
        self.particles = None
        self.time = 0.0
        self.isRunning = False
        self._stopRequested = False

    def setParameters(self, parameters: SimulatorParameters):
        self.parameters = parameters
        self._recreateComponents()
        self.particles = None
        self.time = 0.0

    def _recreateComponents(self):
        self.timeStep = self.parameters.timeStep
        self.gravitationalConstant = self.parameters.gravitationalConstant
        self.calculator = self._createCalculator(self.parameters)
        self.integrator = self._createIntegrator(self.parameters)
        self.particles = self._createDistribution(self.parameters)

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
        if calculatorType == "NEWTON":
            return NewtonCalculator(gravitationalConstant=parameters.gravitationalConstant)
        else:
            return NewtonCalculator(gravitationalConstant=parameters.gravitationalConstant)

    @staticmethod
    def _createDistribution(parameters: SimulatorParameters):
        distributionClasses = {"BASIC": Distribution, "GALAXY": GalaxyDistribution}
        distributionType = parameters.distributionType.upper()
        try:
            return distributionClasses[distributionType]().generate(parameters)
        except KeyError:
            raise ValueError(f"Unknown distribution type: {distributionType}")

    @pyqtSlot()
    def run(self):
        print('NBodySimulator: run() started on thread', QThread.currentThread())
        if QThread.currentThread() != self.thread():
            print("WARNING: run() is NOT running on worker thread!")
        else:
            print("OK: Running on worker thread")
        self.isRunning = True
        self._stopRequested = False
        try:
            if self.particles is None:
                self.particles = self._createDistribution(self.parameters)
            stepCount = 0
            while self.isRunning and not self._stopRequested:
                self.particles = self.integrator.step(self.particles, self.calculator)
                self.time += self.timeStep
                groups = {"main": self.particles}
                self.positionsReady.emit(groups)
                stepCount += 1
        except Exception as e:
            print("Simulation error:", e)
        finally:
            self.isRunning = False
            self.simulationFinished.emit()
            print('NBodySimulator: run() finished')

    @pyqtSlot()
    def stop(self):
        self._stopRequested = True
        self.isRunning = False

    def step(self):
        if self.positions is None:
            self.initialize()
        self.particles = self.integrator.step(self.positions, self.velocities, self.masses, self.calculator)
        self.time += self.timeStep
        return self.positions, self.velocities
