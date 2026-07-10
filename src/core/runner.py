from PyQt5.QtCore import QThread, pyqtSignal
import numpy as cp
from src.core.engine.simulators import NBodySimulator
from src.core.engine.parameters import SimulatorParameters


class SimulationRunner(QThread):
    positionsUpdated = pyqtSignal(dict)   # {group_id: positions_array}
    simulationFinished = pyqtSignal()

    def __init__(self, parameters: SimulatorParameters):
        super().__init__()
        self.parameters = parameters
        self.simulator = None
        self.isRunning = False
        self._stopRequested = False

    def run(self):
        print('Hi')
        self.isRunning = True
        self._stopRequested = False

        try:
            self.simulator = NBodySimulator(self.parameters)
            self.simulator.initialize()
            while self.isRunning and not self._stopRequested:
                positions, _ = self.simulator.step()
                groups = {"main": positions}
                self.positionsUpdated.emit(groups)
                self.msleep(16)

        except Exception as e:
            print("Simulation error:", e)
        finally:
            self.isRunning = False
            self.simulationFinished.emit()

    def stop(self):
        self._stopRequested = True
        self.isRunning = False