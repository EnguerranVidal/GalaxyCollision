import os
import json
import time
from typing import Optional

from PyQt5.QtCore import QUrl, QTimer, Qt, QThread, Q_ARG, QMetaObject, pyqtSlot
from PyQt5.QtGui import QIcon, QDesktopServices
from PyQt5.QtWidgets import *

from core.engine.simulators import NBodySimulator
from src.core.engine.parameters import SimulatorParameters
from src.gui.configEditor import SimulationConfigEditorDock
from src.gui.visualizers.view3d import Universe3dViewWidget
from src.gui.settings import UiSettings, WindowGeometry


class MainWindow(QMainWindow):
    def __init__(self, currentDir: str):
        super().__init__()
        self.icons = {}
        self.settings: Optional[UiSettings] = None
        # FOLDER PATHS & SETTINGS
        self.currentDir = currentDir
        self.settingsPath = os.path.join(self.currentDir, 'settings.json')
        self._checkEnvironment()
        self.loadSettings()

        # SETTING UP USER INTERFACE & RUNNER THREAD
        self.activeParameters = self.settings.parameters
        self.workerThread = QThread(self)
        self.nBodySimulator = NBodySimulator(self.activeParameters)
        self.nBodySimulator.positionsReady.connect(self._onPositionsUpdated)
        self.nBodySimulator.moveToThread(self.workerThread)
        self.workerThread.start()
        self.configEditorDock = SimulationConfigEditorDock(self.activeParameters, self)
        self.configEditorDock.launchSimulationPressed.connect(self._onLaunchSimulation)
        self.configEditorDock.resetSimulationPressed.connect(self._onResetSimulation)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.configEditorDock)
        self.stackedSimulationWidget = QStackedWidget()
        self.simulation3dWidget = Universe3dViewWidget()
        self.stackedSimulationWidget.addWidget(self.simulation3dWidget)
        self.setCentralWidget(self.stackedSimulationWidget)
        print("Main thread:", QThread.currentThread())
        print("Simulator thread:", self.nBodySimulator.thread())

        self._createIcons()
        self._createActions()
        self._createMenuBar()
        self._setupStatusBar()
        self._restoreWindow()

    def _createActions(self):
        # VISIT GITHUB
        self.githubAction = QAction('&Visit GitHub', self)
        self.githubAction.setIcon(self.icons['GITHUB'])
        self.githubAction.setStatusTip('Visit the Project\'s GitHub Repository')
        self.githubAction.triggered.connect(self._openGithub)
        # REPORT ISSUE
        self.reportIssueAction = QAction('&Report Issue', self)
        self.reportIssueAction.setIcon(self.icons['BUG'])
        self.reportIssueAction.setStatusTip('Report an Issue')
        self.reportIssueAction.triggered.connect(self._reportIssue)
        # QUIT APPLICATION
        self.quitAction = QAction('&Quit', self)
        self.quitAction.setStatusTip('Quit Application')
        self.quitAction.triggered.connect(self.close)

    def _createMenuBar(self):
        self.menuBar = self.menuBar()
        ### FILE MENU ###
        self.fileMenu = self.menuBar.addMenu('&File')
        self.fileMenu.addAction(self.quitAction)
        ### HELP MENU ###
        self.helpMenu = self.menuBar.addMenu('&Help')
        self.helpMenu.addAction(self.githubAction)
        self.helpMenu.addAction(self.reportIssueAction)

    def _createIcons(self):
        self.iconPath = os.path.join(self.currentDir, f'src/assets/icons')
        self.icons['BUG'] = QIcon(os.path.join(self.iconPath, 'bug.png'))
        self.icons['GITHUB'] = QIcon(os.path.join(self.iconPath, 'github.png'))

    @pyqtSlot(SimulatorParameters)
    def _onLaunchSimulation(self, parameters: SimulatorParameters):
        self.activeParameters = parameters
        QMetaObject.invokeMethod(self.nBodySimulator, "stop", Qt.QueuedConnection)
        QTimer.singleShot(50, lambda: self._startSimulation(parameters))

    @pyqtSlot()
    def _onResetSimulation(self):
        QMetaObject.invokeMethod(self.nBodySimulator, "stop", Qt.QueuedConnection)
        QTimer.singleShot(50, lambda: self._startSimulation(self.activeParameters))

    def _startSimulation(self, parameters: SimulatorParameters):
        self.nBodySimulator.setParameters(parameters)
        QMetaObject.invokeMethod(self.nBodySimulator, "run", Qt.QueuedConnection)

    def _onSimulationFinished(self):
        self.statusBar().showMessage("Simulation finished", 3000)

    def _onPositionsUpdated(self, groups: dict):
        if self.simulation3dWidget:
            # self.simulation3dWidget.updateData(groups)
            self._updateFps()

    def _checkEnvironment(self):
        if not os.path.exists(self.settingsPath):
            settings = UiSettings().toDict()
            with open(self.settingsPath, 'w') as f:
                json.dump(settings, f)

    def loadSettings(self):
        with open(self.settingsPath, 'r') as f:
            settings = json.load(f)
        self.settings = UiSettings.fromDict(settings)

    def saveSettings(self):
        if self.settings is None:
            return
        settings = self.settings.toDict()
        with open(self.settingsPath, "w") as f:
            json.dump(settings, f)

    def _restoreWindow(self):
        self.setWindowTitle('Galaxy Collision')
        if self.settings.window.maximized:
            self.showMaximized()
        else:
            windowGeometry = self.settings.window.geometry
            self.setGeometry(windowGeometry.x, windowGeometry.y, windowGeometry.width, windowGeometry.height)

    def _center(self):
        frameGeometry = self.frameGeometry()
        screenCenter = QDesktopWidget().availableGeometry().center()
        frameGeometry.moveCenter(screenCenter)
        self.move(frameGeometry.topLeft())

    @staticmethod
    def _openGithub():
        QDesktopServices.openUrl(QUrl("https://github.com/EnguerranVidal/GalaxyCollision"))

    @staticmethod
    def _reportIssue():
        QDesktopServices.openUrl(QUrl("https://github.com/EnguerranVidal/GalaxyCollision/issues/new"))

    def _setupStatusBar(self):
        self.lastUpdate = time.perf_counter()
        self.avgFps = 0.0
        self.fpsLabel = QLabel('Fps : ---')
        self.fpsLabel.setStyleSheet('border: 0;')
        self.statusBar().addPermanentWidget(self.fpsLabel)
        self.statusBar().showMessage('Ready')
        self.statusDateTimer = QTimer()
        self.statusDateTimer.timeout.connect(self._updateStatus)
        self.statusDateTimer.start(1000)

    def _updateFps(self):
        now = time.perf_counter()
        fps = 1.0 / (now - self.lastUpdate)
        self.lastUpdate = now
        self.avgFps = self.avgFps * 0.8 + fps * 0.2
        self.fpsLabel.setText(f'FPS : {self.avgFps:.1f}')

    def _updateStatus(self):
        self.statusBar().showMessage(f'Time: {getattr(self.nBodySimulator, "time", 0):.1f}s')

    def closeEvent(self, event):
        if self.nBodySimulator and self.nBodySimulator.isRunning:
            self.nBodySimulator.stop()
            if not self.nBodySimulator.wait(1000):
                print("Warning: Simulator did not stop gracefully")
        self.settings.window.maximized = self.isMaximized()
        if not self.isMaximized():
            g = self.geometry()
            self.settings.window.geometry = WindowGeometry.fromDict({'X': g.x(), 'Y': g.y(), 'WIDTH': g.width(), 'HEIGHT': g.height()})
        self.saveSettings()
        self.workerThread.quit()
        self.workerThread.wait(1000)
        event.accept()