from __future__ import annotations
import time
import numpy as np
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

import engine
from src.gui.parameters import *


class NBodySimulator(QObject):
    positionsReady = pyqtSignal(dict)
    simulationFinished = pyqtSignal()

    MAX_FPS = 30.0

    def __init__(self, parameters, parent=None):
        super().__init__(parent)
        self.parameters = parameters
        self.cppParameters = None
        self.particles = None
        self.distribution = None
        self.calculator = None
        self.integrator = None
        self.isRunning = False
        self.simulationTime = 0.0
        self.frameInterval = 1.0 / self.MAX_FPS

    def initialize(self):
        self.cppParameters = self.parameters.toCpp()
        self.distribution = self._createDistribution(self.cppParameters)
        self.calculator = self._createCalculator(self.cppParameters)
        self.integrator = self._createIntegrator(self.cppParameters)
        self.particles = self.distribution.generate(self.cppParameters)
        self.simulationTime = 0.0

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

    @pyqtSlot(SimulatorParameters)
    def prepareAndRun(self, parameters):
        self.parameters = parameters
        self.particles = None
        self.initialize()
        self.run()

    @pyqtSlot()
    def run(self):
        if self.particles is None:
            self.initialize()
        self.isRunning = True
        nextFrame = time.perf_counter()
        while self.isRunning:
            self.integrator.step(self.particles, self.calculator)
            self.simulationTime += self.parameters.timeStep
            now = time.perf_counter()
            if now >= nextFrame:
                self.positionsReady.emit({"default": self._positionsAsArray()})
                nextFrame = now + self.frameInterval
            if not self.cppParameters.endless and self.simulationTime >= self.cppParameters.maxTime:
                break
            sleepFor = nextFrame - time.perf_counter()
            if sleepFor > 0:
                time.sleep(sleepFor)
        self.isRunning = False
        self.simulationFinished.emit()

    @pyqtSlot()
    def step(self):
        if self.particles is None:
            self.initialize()
        self.integrator.step(self.particles, self.calculator)
        self.simulationTime += self.cppParameters.timeStep
        self.positionsReady.emit({"default": self._positionsAsArray()})

    @pyqtSlot()
    def stop(self):
        self.isRunning = False

    def _positionsAsArray(self):
        positions = self.particles.getPositions()
        if not positions:
            return np.zeros((0, 3), dtype=np.float32)
        return np.asarray([[position.x, position.y, position.z] for position in positions], dtype=np.float32)
