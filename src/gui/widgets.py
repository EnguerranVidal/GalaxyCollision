import os

from PyQt5.QtCore import Qt, pyqtSignal, QSize
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import *


class SimulationTimeBar(QWidget):
    pauseToggled = pyqtSignal(bool)
    speedChanged = pyqtSignal(float)

    def __init__(self, parent=None, currentDirectory:str = ''):
        super().__init__(parent)
        self.iconPath = os.path.join(currentDirectory, 'src/assets/icons')
        self.allowedSpeeds = [0.1, 0.25, 0.5, 1, 2, 4, 8, 16]
        self.speedIndex = self.allowedSpeeds.index(1)
        self._paused = False

        # WIDGETS & LAYOUT
        self.slowDownButton = SquareIconButton(os.path.join(self.iconPath, 'slow-down.png'))
        self.slowDownButton.setToolTip("Slow down simulation")
        self.slowDownButton.clicked.connect(self._slow)
        self.playPauseButton = SquareIconButton(os.path.join(self.iconPath, 'pause.png'))
        self.playPauseButton.setToolTip("Pause / resume simulation")
        self.playPauseButton.setCheckable(True)
        self.playPauseButton.toggled.connect(self._onPauseToggled)
        self.fastForwardButton = SquareIconButton(os.path.join(self.iconPath, 'fast-forward.png'))
        self.fastForwardButton.setToolTip("Fast forward simulation")
        self.fastForwardButton.clicked.connect(self._fast)
        self.speedButton = QPushButton("x1")
        self.speedButton.setMinimumWidth(40)
        self.speedButton.clicked.connect(self._resetSpeed)
        self.simulationTimeLabel = QLabel("Time: 0.000")
        self.simulationTimeLabel.setMinimumWidth(100)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        layout.addWidget(self.simulationTimeLabel)
        layout.addStretch(1)
        layout.addWidget(self.slowDownButton)
        layout.addWidget(self.playPauseButton)
        layout.addWidget(self.fastForwardButton)
        layout.addWidget(self.speedButton)
        layout.addStretch(1)
        self.setEnabled(False)

    def setPaused(self, paused: bool):
        self._paused = bool(paused)
        self.playPauseButton.blockSignals(True)
        self.playPauseButton.setChecked(self._paused)
        icon = "play.png" if self._paused else "pause.png"
        self.playPauseButton.setIconPath(os.path.join(self.iconPath, icon))
        self.playPauseButton.setToolTip("Resume simulation" if self._paused else "Pause simulation")
        self.playPauseButton.blockSignals(False)

    def setSpeed(self, speed: float):
        speed = float(speed)
        self.speedIndex = min(range(len(self.allowedSpeeds)), key=lambda i: abs(self.allowedSpeeds[i] - speed))
        self.speedButton.setText(f"×{self.allowedSpeeds[self.speedIndex]:g}")

    def _onPauseToggled(self, paused: bool):
        self._paused = paused
        icon = "play.png" if paused else "pause.png"
        self.playPauseButton.setIconPath(os.path.join(self.iconPath, icon))
        self.playPauseButton.setToolTip("Resume" if paused else "Pause")
        self.pauseToggled.emit(paused)

    def _slow(self):
        if self.speedIndex > 0:
            self.speedIndex -= 1
            speed = self.allowedSpeeds[self.speedIndex]
            self.speedChanged.emit(speed)
            self.speedButton.setText(f"x{speed:g}")

    def _fast(self):
        if self.speedIndex < len(self.allowedSpeeds) - 1:
            self.speedIndex += 1
            speed = self.allowedSpeeds[self.speedIndex]
            self.speedChanged.emit(speed)
            self.speedButton.setText(f"x{speed:g}")

    def _resetSpeed(self):
        self.speedIndex = self.allowedSpeeds.index(1)
        speed = self.allowedSpeeds[self.speedIndex]
        self.speedChanged.emit(speed)
        self.speedButton.setText(f"x{speed:g}")

    def setSimulationTime(self, simulationTime: float):
        self.simulationTimeLabel.setText(f"Time: {simulationTime:.3f}")



class SquareIconButton(QPushButton):
    def __init__(self, iconPath: str, parent=None, size=24, flat=False):
        super(SquareIconButton, self).__init__(parent)
        self.iconPath = iconPath
        self.setIcon(QIcon(self.iconPath))
        self.setIconSize(QSize(size, size))
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        if flat:
            styleSheet = 'border: none; text-align: center;'
            self.setFlat(True)
        else:
            styleSheet = 'text-align: center;'
        self.setStyleSheet(styleSheet)
        self.setAutoFillBackground(False)

    def setIconPath(self, iconPath: str):
        self.iconPath = iconPath
        self.setIcon(QIcon(iconPath))

    def setIconObject(self, icon: QIcon):
        self.iconPath = None
        self.setIcon(icon)

    def setIconSize(self, size):
        super().setIconSize(size)
        self.setFixedSize(size)

    def sizeHint(self):
        return self.iconSize()