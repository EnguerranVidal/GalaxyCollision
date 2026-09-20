from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict
import time
import numpy as np
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

import engine
from src.gui.solver.parameters import SolverParameters
from src.gui.solver.distributions import generate


class NBodySolver(QObject):
    positionsReady = pyqtSignal(object)
    simulationFinished = pyqtSignal()

    MAX_FPS = 60.0

    def __init__(self, parameters, parent=None):
        super().__init__(parent)
        self.parameters = parameters
        self.cppParameters = None
        self.particles = None
        self.calculator = None
        self.integrator = None
        self.isRunning = False
        self.isPaused = False
        self.simulationTime = 0.0

    def initialize(self):
        self.cppParameters = self.parameters.toCpp()
        self.calculator = self._createCalculator(self.cppParameters)
        self.integrator = self._createIntegrator(self.cppParameters)
        self.particles = generate(self.parameters)
        self.simulationTime = 0.0

    @staticmethod
    def _createCalculator(parameters: engine.SolverParameters):
        calculatorType = parameters.calculatorType.upper()
        if calculatorType == "BARNES_HUT":
            calculator = engine.BarnesHutCalculator(parameters.theta, parameters.gravitationalConstant)
        else:
            calculator = engine.NewtonCalculator(parameters.gravitationalConstant)
        return calculator

    @staticmethod
    def _createIntegrator(parameters: engine.SolverParameters):
        integratorType = parameters.integratorType.upper()
        if integratorType == "RK4":
            integrator = engine.RK4Integrator(parameters.timeStep)
        else:
            integrator = engine.EulerIntegrator(parameters.timeStep)
        return integrator

    @pyqtSlot(SolverParameters)
    def prepareAndRun(self, parameters):
        self.parameters = parameters
        self.particles = None
        self.initialize()
        self.run()

    def _emitState(self):
        state = State(
            time=self.simulationTime,
            positions={"default": self._positionsAsArray()},
            velocities={"default": self._velocitiesAsArray()},
            accelerations={"default": self._accelerationsAsArray()},
            masses={"default": self._massesAsArray()},
            massCenter=self._massCenterAsArray(),
        )
        self.positionsReady.emit(state)

    @pyqtSlot()
    def run(self):
        if self.particles is None:
            self.initialize()
        self.isRunning, self.isPaused = True, False
        self.simulationTime = 0.0
        displayTimeStep = self._displayTimeStep()
        nextDisplayTime = displayTimeStep
        self._emitState()
        while self.isRunning:
            if self.isPaused:
                time.sleep(0.01)
                continue
            wallFrameStart = time.perf_counter()
            while self.simulationTime < nextDisplayTime and self.isRunning and not self.isPaused:
                calculationTimeStep = self._calculationTimeStep()
                remainingTime = nextDisplayTime - self.simulationTime
                if calculationTimeStep > remainingTime:
                    calculationTimeStep = remainingTime
                if calculationTimeStep <= 0.0:
                    break
                self.integrator.step(self.particles, self.calculator)
                self.simulationTime += float(self.parameters.timeStep)
            if not self.isRunning:
                break
            self._emitState()
            nextDisplayTime += self._displayTimeStep()
            if not self.cppParameters.endless and self.simulationTime >= self.cppParameters.maxTime:
                break
            frameWall = 1.0 / max(float(self.parameters.maxFps), 1e-6)
            sleepFor = frameWall - (time.perf_counter() - wallFrameStart)
            if sleepFor > 0.0:
                time.sleep(sleepFor)
        self.isRunning = False
        self.isPaused = False
        self.simulationFinished.emit()

    @pyqtSlot()
    def step(self):
        if self.particles is None:
            self.initialize()
        self.integrator.step(self.particles, self.calculator)
        self.simulationTime += self.cppParameters.timeStep
        state = State(time=self.simulationTime,
                      positions={"default": self._positionsAsArray()},
                      velocities={"default": self._velocitiesAsArray()},
                      accelerations={"default": self._accelerationsAsArray()},
                      masses={"default": self._massesAsArray()},
                      massCenter=self._massCenterAsArray())
        self.positionsReady.emit(state)

    @pyqtSlot()
    def pause(self):
        self.isPaused = True

    @pyqtSlot()
    def resume(self):
        self.isPaused = False

    @pyqtSlot()
    def stop(self):
        self.isPaused = False
        self.isRunning = False

    def _positionsAsArray(self):
        n = self.particles.getNbParticles()
        outputArray = np.empty((n, 3), dtype=np.float32)
        self.particles.copyPositionsTo(outputArray)
        return outputArray

    def _velocitiesAsArray(self):
        n = self.particles.getNbParticles()
        outputArray = np.empty((n, 3), dtype=np.float32)
        self.particles.copyVelocitiesTo(outputArray)
        return outputArray

    def _accelerationsAsArray(self):
        n = self.particles.getNbParticles()
        outputArray = np.empty((n, 3), dtype=np.float32)
        self.particles.copyAccelerationsTo(outputArray)
        return outputArray

    def _massesAsArray(self):
        n = self.particles.getNbParticles()
        outputArray = np.empty(n, dtype=np.float32)
        self.particles.copyMassesTo(outputArray)
        return outputArray

    def _massCenterAsArray(self):
        massCenter, _ = self.particles.massCenter()
        if not massCenter:
            return np.zeros(3, dtype=np.float32)
        return np.asarray([massCenter.x, massCenter.y, massCenter.z], dtype=np.float32)

    def _displayTimeStep(self):
        maxFps = max(float(self.parameters.maxFps), 1e-6)
        simulationSpeed = max(float(self.parameters.simulationSpeed), 0.0)
        return simulationSpeed / maxFps

    def _calculationTimeStep(self):
        return self.parameters.timeStep

    @pyqtSlot(float)
    def setSimulationSpeed(self, speed: float):
        self.parameters.simulationSpeed = max(float(speed), 0.0)


@dataclass
class State:
    time: float = 0.0
    positions: Dict[str, np.ndarray] = field(default_factory=dict)
    velocities: Dict[str, np.ndarray] = field(default_factory=dict)
    accelerations: Dict[str, np.ndarray] = field(default_factory=dict)
    masses: Dict[str, np.ndarray] = field(default_factory=dict)
    massCenter: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float32))
