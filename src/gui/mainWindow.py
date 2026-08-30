import os
import json
import time

from PyQt5.QtCore import QUrl, QTimer, Qt, QThread, Q_ARG, QMetaObject, pyqtSlot, QEventLoop
from PyQt5.QtGui import QIcon, QDesktopServices
from PyQt5.QtWidgets import *

from src.gui.solver.simulator import NBodySimulator, State
from src.gui.solver.parameters import SimulatorParameters
from src.gui.configEditor import SimulationConfigEditorDock
from src.gui.visualizers.view3d import Universe3dViewWidget
from src.gui.settings import UiSettings


class MainWindow(QMainWindow):
    def __init__(self, currentDir: str):
        super().__init__()
        self.icons = {}
        self.settings = UiSettings()
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
        self.nBodySimulator.simulationFinished.connect(self._onSimulationFinished)
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
        self.simulation3dWidget.setViewSettings(self.settings.view)
        print("Main thread:", QThread.currentThread())
        print("Simulator thread:", self.nBodySimulator.thread())

        self._createIcons()
        self._createActions()
        self._createMenuBar()
        self._setupStatusBar()
        self._restoreWindow()

    def _createActions(self):
        # SHOW BARYCENTER
        self.showBarycenterAction = QAction('&Barycenter', self)
        self.showBarycenterAction.setIcon(self.icons['CENTER_GRAVITY'])
        self.showBarycenterAction.setStatusTip('Show the Simulation\'s Barycenter.')
        self.showBarycenterAction.setCheckable(True)
        self.showBarycenterAction.setChecked(self.settings.view.showBarycenter)
        self.showBarycenterAction.toggled.connect(self._toggleBarycenter)
        self.showBarycenterAction.setIconVisibleInMenu(False)
        # CENTER ON BARYCENTER
        self.centerOnBarycenterAction = QAction('&Center on Barycentre', self)
        self.centerOnBarycenterAction.setIcon(self.icons['ARROWS_CENTER'])
        self.centerOnBarycenterAction.setStatusTip('Center Camera on Barycenter.')
        self.centerOnBarycenterAction.setCheckable(True)
        self.centerOnBarycenterAction.setChecked(self.settings.view.centerOnBarycenter)
        self.centerOnBarycenterAction.toggled.connect(self._toggleCenterOnBarycenter)
        self.centerOnBarycenterAction.setIconVisibleInMenu(False)
        # SHOW VELOCITY VECTORS
        self.showVelocityVectorsAction = QAction("&Velocity vectors", self)
        self.showVelocityVectorsAction.setStatusTip('Show Velocity Vectors.')
        self.showVelocityVectorsAction.setCheckable(True)
        self.showVelocityVectorsAction.setChecked(self.settings.view.showVelocityVectors)
        self.showVelocityVectorsAction.toggled.connect(self._toggleVelocityVectors)
        self.showVelocityVectorsAction.setIconVisibleInMenu(False)
        # SHOW ACCELERATION VECTORS
        self.showAccelerationVectorsAction = QAction("&Acceleration vectors", self)
        self.showAccelerationVectorsAction.setStatusTip('Show Acceleration Vectors.')
        self.showAccelerationVectorsAction.setCheckable(True)
        self.showAccelerationVectorsAction.setChecked(self.settings.view.showAccelerationVectors)
        self.showAccelerationVectorsAction.toggled.connect(self._toggleAccelerationVectors)
        self.showAccelerationVectorsAction.setIconVisibleInMenu(False)
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
        ### VIEW MENU ###
        self.viewMenu = self.menuBar.addMenu('&View')
        self.viewMenu.addAction(self.showBarycenterAction)
        self.viewMenu.addAction(self.centerOnBarycenterAction)
        self.viewMenu.addSeparator()
        self.viewMenu.addAction(self.showVelocityVectorsAction)
        self.viewMenu.addAction(self.showAccelerationVectorsAction)
        ### HELP MENU ###
        self.helpMenu = self.menuBar.addMenu('&Help')
        self.helpMenu.addAction(self.githubAction)
        self.helpMenu.addAction(self.reportIssueAction)

    def _createIcons(self):
        self.iconPath = os.path.join(self.currentDir, f'src/assets/icons')
        self.icons['ARROWS_CENTER'] = QIcon(os.path.join(self.iconPath, 'arrows-center.png'))
        self.icons['BUG'] = QIcon(os.path.join(self.iconPath, 'bug.png'))
        self.icons['CENTER_GRAVITY'] = QIcon(os.path.join(self.iconPath, 'center-gravity.png'))
        self.icons['GITHUB'] = QIcon(os.path.join(self.iconPath, 'github.png'))

    @pyqtSlot(SimulatorParameters)
    def _onLaunchSimulation(self, parameters: SimulatorParameters):
        self.activeParameters = parameters
        self.settings.parameters = parameters
        self.saveSettings()
        self._restartSimulation(parameters)

    @pyqtSlot()
    def _onResetSimulation(self):
        self._restartSimulation(self.activeParameters)

    def _restartSimulation(self, parameters: SimulatorParameters):
        self.nBodySimulator.stop()
        QTimer.singleShot(30, lambda: self._startNewSimulation(parameters))

    def _startNewSimulation(self, parameters: SimulatorParameters):
        QMetaObject.invokeMethod(self.nBodySimulator, "prepareAndRun", Qt.QueuedConnection, Q_ARG(object, parameters))

    def _onSimulationFinished(self):
        self.statusBar().showMessage("Simulation finished", 3000)

    def _onPositionsUpdated(self, state: State):
        if self.simulation3dWidget:
            self.simulation3dWidget.updateState(state)
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
            windowSettings = self.settings.window
            self.setGeometry(windowSettings.x, windowSettings.y, windowSettings.width, windowSettings.height)

    def _center(self):
        frameGeometry = self.frameGeometry()
        screenCenter = QDesktopWidget().availableGeometry().center()
        frameGeometry.moveCenter(screenCenter)
        self.move(frameGeometry.topLeft())

    def _toggleBarycenter(self, checked: bool):
        self.settings.view.showBarycenter = checked
        self.simulation3dWidget.setShowBarycenter(checked)
        self.saveSettings()

    def _toggleCenterOnBarycenter(self, checked: bool):
        self.settings.view.centerOnBarycenter = checked
        self.simulation3dWidget.setCenterOnBarycenter(checked)
        self.saveSettings()

    def _toggleVelocityVectors(self, checked: bool):
        self.settings.view.showVelocityVectors = checked
        self.simulation3dWidget.velocityVectorLength = self.settings.view.velocityVectorLength
        self.simulation3dWidget.referenceVelocity = self.settings.view.referenceVelocity
        self.simulation3dWidget.setShowVelocityVectors(checked)
        self.saveSettings()

    def _toggleAccelerationVectors(self, checked: bool):
        self.settings.view.showAccelerationVectors = checked
        self.simulation3dWidget.accelerationVectorLength = self.settings.view.accelerationVectorLength
        self.simulation3dWidget.referenceAcceleration = self.settings.view.referenceAcceleration
        self.simulation3dWidget.setShowAccelerationVectors(checked)
        self.saveSettings()

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

    def _updateFps(self):
        now = time.perf_counter()
        fps = 1.0 / (now - self.lastUpdate)
        self.lastUpdate = now
        self.avgFps = self.avgFps * 0.8 + fps * 0.2
        self.fpsLabel.setText(f'FPS : {self.avgFps:.1f}')

    def closeEvent(self, event):
        if self.nBodySimulator and self.nBodySimulator.isRunning:
            print("Closing: Stopping simulation...")
            self.nBodySimulator.stop()
            loop = QEventLoop()
            QTimer.singleShot(300, loop.quit)
            finished = False

            def onFinished():
                nonlocal finished
                finished = True
                loop.quit()

            connection = self.nBodySimulator.simulationFinished.connect(onFinished)
            loop.exec_()
            self.nBodySimulator.simulationFinished.disconnect(connection)
            if not finished:
                print("Warning: Simulator did not stop gracefully within timeout")
        self.settings.parameters = self.configEditorDock.getUiParameters()
        self.settings.window.maximized = self.isMaximized()
        if not self.isMaximized():
            g = self.geometry()
            self.settings.window.x, self.settings.window.y = g.x(), g.y()
            self.settings.window.width, self.settings.window.height = g.width(), g.height()
        self.saveSettings()
        if self.workerThread and self.workerThread.isRunning():
            print("Closing: Quitting worker thread...")
            self.workerThread.quit()
            if not self.workerThread.wait(1500):
                print("Warning: Worker thread did not quit gracefully - forcing termination")
                self.workerThread.terminate()
                self.workerThread.wait(500)
        event.accept()
