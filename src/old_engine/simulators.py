import time

from PyQt5.QtCore import QThread, pyqtSignal, QObject, pyqtSlot

from src.engine.integrators import EulerExplicit, RK4
from src.engine.calculators import NewtonCalculator, BarnesHutCalculator
from src.engine.distributions import Distribution, GalaxyDistribution
from src.engine.parameters import SimulatorParameters


class NBodySimulator(QObject):
    positionsReady = pyqtSignal(dict)
    simulationFinished = pyqtSignal()

    def __init__(self, parameters: SimulatorParameters):
        super().__init__()
        self.particles = None
        self.parameters = parameters
        self._initializeComponents()
        self.time = 0.0
        self.isRunning = False
        self._stopRequested = False

    def _initializeComponents(self):
        self.timeStep = self.parameters.timeStep
        self.gravitationalConstant = self.parameters.gravitationalConstant
        self.calculator = self._createCalculator(self.parameters)
        self.integrator = self._createIntegrator(self.parameters)
        self.particles = self._createDistribution(self.parameters)

    def setParameters(self, parameters: SimulatorParameters):
        self.parameters = parameters
        self._initializeComponents()
        self.time = 0.0

    @staticmethod
    def _createIntegrator(parameters: SimulatorParameters):
        integratorType = parameters.integratorType.upper().replace("_", "").replace("-", "")
        if integratorType == "EULER_EXPLICIT":
            return EulerExplicit(parameters.timeStep)
        else:
            return RK4(parameters.timeStep)

    @staticmethod
    def _createCalculator(parameters: SimulatorParameters):
        calculatorType = parameters.calculatorType.upper().replace("_", "").replace("-", "")
        if calculatorType == 'NEWTON':
            return NewtonCalculator(gravitationalConstant=parameters.gravitationalConstant)
        elif calculatorType == 'BARNES_HUT':
            return BarnesHutCalculator(gravitationalConstant=parameters.gravitationalConstant)
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

    @pyqtSlot(SimulatorParameters)
    def prepareAndRun(self, parameters: SimulatorParameters):
        print('NBodySimulator: prepareAndRun on thread', QThread.currentThread())
        self.parameters = parameters
        self._initializeComponents()
        self.isRunning = True
        self._stopRequested = False
        self._runLoop()

    @pyqtSlot()
    def _runLoop(self):
        print('NBodySimulator: _runLoop() started on thread', QThread.currentThread())
        stepCount = 0
        lastEmitTime = 0.0
        minimumEmitInterval = 1.0 / 60.0
        try:
            while self.isRunning and not self._stopRequested:
                self.step()
                now = time.perf_counter()
                if now - lastEmitTime >= minimumEmitInterval:
                    groups = {"main": self._positionsForDisplay()}
                    if groups["main"] is not None:
                        self.positionsReady.emit(groups)
                    lastEmitTime = now
                stepCount += 1
                if not self.parameters.endless and self.time >= self.parameters.maxTime:
                    break
        except Exception as e:
            print("Simulation error:", e)
        finally:
            self.isRunning = False
            self.simulationFinished.emit()
            print(f'NBodySimulator: _runLoop finished after {stepCount} steps')

    @pyqtSlot()
    def stop(self):
        print("NBodySimulator: stop() called")
        self._stopRequested = True
        self.isRunning = False

    def step(self):
        self.particles = self.integrator.step(self.particles, self.calculator)
        self.time += self.timeStep

    def _positionsForDisplay(self):
        positions = self.particles.positions
        if self.particles.device == "gpu":
            return positions.get()
        return positions
