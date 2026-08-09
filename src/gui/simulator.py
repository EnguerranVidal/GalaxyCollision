from __future__ import annotations
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

import engine
from src.gui.parameters import *


class NBodySimulator(QObject):
    positionsReady = pyqtSignal(dict)
    simulationFinished = pyqtSignal()

    def __init__(self, parameters, parent=None):
        super().__init__(parent)
        print(engine.__file__)
        self.parameters = parameters
        self.cppParameters = None
        self.particles = None
        self.distribution = None
        self.calculator = None
        self.integrator = None
        self._running = False
        self._paused = False
        self._time = 0.0

    def initialize(self):
        self.cppParameters = self.parameters.toCpp()
        self.distribution = self._createDistribution(self.cppParameters)
        self.calculator = self._createCalculator(self.cppParameters)
        self.integrator = self._createIntegrator(self.cppParameters)
        self.particles = self.distribution.generate(self.cppParameters)
        self._time = 0.0

    @staticmethod
    def _createDistribution(parameters: engine.SimulatorParameters):
        distributionType = parameters.distributionType.upper()
        if distributionType == "GALAXY":
            distribution = engine.GalaxyDistribution()
        else:
            distribution = engine.BasicDistribution()
        return distribution


    @staticmethod
    def _createCalculator(parameters: engine.SimulatorParameters):
        calculatorType = parameters.calculatorType.upper()
        if calculatorType == "BARNES_HUT":
            calculator = engine.BarnesHutCalculator(parameters.theta, parameters.gravitationalConstant)
        else:
            calculator = engine.NewtonCalculator(parameters.gravitationalConstant)
        return calculator

    @staticmethod
    def _createIntegrator(parameters: engine.SimulatorParameters):
        integratorType = parameters.integratorType.upper()
        if integratorType == "RK4":
            integrator = engine.RK4Integrator(parameters.timeStep)
        else:
            integrator = engine.EulerIntegrator(parameters.timeStep)
        return integrator

    @pyqtSlot()
    def run(self):
        if self.particles is None:
            self.initialize()
        self._running = True
        while self._running:
            if self._paused:
                continue
            self.integrator.step(self.particles, self.calculator)
            self._time += self.parameters.timeStep
            self.positionsReady.emit({"time": self._time, "positions": self.particles.getPositions()})
            if not self.cppParameters.endless and self._time >= self.cppParameters.maxTime:
                break
        self._running = False
        self.simulationFinished.emit()

    @pyqtSlot()
    def step(self):
        if self.particles is None:
            self.initialize()
        self.integrator.step(self.particles, self.calculator)
        self._time += self.cppParameters.timeStep
        self.positionsReady.emit({"time": self._time, "positions": self.particles.getPositions()})
