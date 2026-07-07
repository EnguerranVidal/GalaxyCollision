from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import *

from src.core.engine.parameters import SimulatorParameters, BasicDistributionParameters, GalaxyDistributionParameters


class SimulationConfigDock(QDockWidget):
    parametersChanged = pyqtSignal(SimulatorParameters)

    def __init__(self, initialParameters: SimulatorParameters = None, parent=None):
        super().__init__("Simulation Configuration", parent)
        self.parameters = initialParameters or SimulatorParameters()
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.setFloating(False)

        # MAIN SIMULATION PARAMETERS
        self.simulationTab = QWidget()
        simulationInfoGroup = QGroupBox("Simulation Information")
        self.nameEdit = QLineEdit("My N-Body Simulation")
        self.saveResultsCheck = QCheckBox("Save Results to Disk")
        simulationInfoForm = QFormLayout(simulationInfoGroup)
        simulationInfoForm.addRow("Name:", self.nameEdit)
        simulationInfoForm.addRow("", self.saveResultsCheck)


        simulationLayout = QVBoxLayout(self.simulationTab)
        simulationLayout.addWidget(simulationInfoGroup)

        # DISTRIBUTION PARAMETERS
        self.distributionTab = QStackedWidget()
        self.basicWidget = BasicDistributionWidget()
        self.galaxyWidget = GalaxyDistributionWidget()
        self.distributionTab.addWidget(self.basicWidget)
        self.distributionTab.addWidget(self.galaxyWidget)

        # MAIN LAYOUT & CONTAINER
        self.container = QWidget()
        self.editorTabs = QTabWidget()
        self.editorTabs.addTab(self.simulationTab, "Simulation")
        self.editorTabs.addTab(self.distributionTab, "Distribution")
        self.reinitializationButton = QPushButton("Reinitialize Simulation")
        self.reinitializationButton.clicked.connect(self._reinitialize)
        mainLayout = QVBoxLayout(self.container)
        mainLayout.addWidget(self.editorTabs)
        mainLayout.addWidget(self.reinitializationButton)
        self.setWidget(self.container)


class BasicDistributionWidget(QWidget):
    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or BasicDistributionParameters()


class GalaxyDistributionWidget(QWidget):
    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or GalaxyDistributionParameters()