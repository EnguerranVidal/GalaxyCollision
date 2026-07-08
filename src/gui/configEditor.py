from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import *

from src.core.engine.parameters import SimulatorParameters, BasicDistributionParameters, GalaxyDistributionParameters


class SimulationConfigDock(QDockWidget):
    launchSimulation = pyqtSignal(SimulatorParameters)

    def __init__(self, initialParameters: SimulatorParameters = None, parent=None):
        super().__init__("Simulation Configuration", parent)
        self.activeParameters = initialParameters or SimulatorParameters()
        self.uiParameters = self.activeParameters
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.setFloating(False)

        # MAIN SIMULATION PARAMETERS
        self.simulationTab = QWidget()
        self.nameEdit = QLineEdit("My N-Body Simulation")
        self.saveResultsCheck = QCheckBox("Save Results to Disk")
        simulationInfoGroup = QGroupBox("Simulation Information")
        simulationInfoForm = QFormLayout(simulationInfoGroup)
        simulationInfoForm.addRow("Name:", self.nameEdit)
        simulationInfoForm.addRow("", self.saveResultsCheck)
        self.distributionTypeCombo = QComboBox()
        self.distributionTypeCombo.addItems(["Galaxy", "Basic"])
        self.distributionTypeCombo.currentTextChanged.connect(self._onDistributionTypeChanged)
        self.calculatorTypeCombo = QComboBox()
        self.calculatorTypeCombo.addItems(["Barnes-Hut", "Newton"])
        self.calculatorTypeCombo.currentTextChanged.connect(self._onUIChanged)
        self.integratorTypeCombo = QComboBox()
        self.integratorTypeCombo.addItems(["RK4", "EulerExplicit"])
        self.integratorTypeCombo.currentTextChanged.connect(self._onUIChanged)
        self.gravitationalConstantSpin = QDoubleSpinBox()
        self.gravitationalConstantSpin.setRange(0.01, 100.0)
        self.gravitationalConstantSpin.setDecimals(4)
        self.gravitationalConstantSpin.valueChanged.connect(self._onUIChanged)
        self.timeStepSpin = QDoubleSpinBox()
        self.timeStepSpin.setRange(0.00001, 0.5)
        self.timeStepSpin.setDecimals(6)
        self.timeStepSpin.valueChanged.connect(self._onUIChanged)
        self.thetaValueSpin = QDoubleSpinBox()
        self.thetaValueSpin.setRange(0.1, 1.0)
        self.thetaValueSpin.setSingleStep(0.05)
        self.thetaValueSpin.valueChanged.connect(self._onUIChanged)
        simulationCoreGroup = QGroupBox("Core Simulation Parameters")
        simulationCoreForm = QFormLayout(simulationCoreGroup)
        simulationCoreForm.addRow("Distribution Type:", self.distributionTypeCombo)
        simulationCoreForm.addRow("Force Calculator:", self.calculatorTypeCombo)
        simulationCoreForm.addRow("Integrator:", self.integratorTypeCombo)
        simulationCoreForm.addRow("Gravitational Constant (G):", self.gravitationalConstantSpin)
        simulationCoreForm.addRow("Time Step (dt):", self.timeStepSpin)
        simulationCoreForm.addRow("Barnes-Hut Theta:", self.thetaValueSpin)
        self.endlessCheck = QCheckBox("Run Indefinitely")
        self.endlessCheck.stateChanged.connect(self._onUIChanged)
        self.maxTimeSpin = QDoubleSpinBox()
        self.maxTimeSpin.setRange(10, 100000)
        self.maxTimeSpin.valueChanged.connect(self._onUIChanged)
        simulationRuntimeGroup = QGroupBox("Runtime")
        simulationRuntimeForm = QFormLayout(simulationRuntimeGroup)
        simulationRuntimeForm.addRow("", self.endlessCheck)
        simulationRuntimeForm.addRow("Max Time:", self.maxTimeSpin)
        simulationLayout = QVBoxLayout(self.simulationTab)
        simulationLayout.addWidget(simulationInfoGroup)
        simulationLayout.addWidget(simulationCoreGroup)
        simulationLayout.addWidget(simulationRuntimeGroup)

        # DISTRIBUTION PARAMETERS
        self.distributionTab = QStackedWidget()
        self.basicWidget = BasicDistributionWidget()
        self.galaxyWidget = GalaxyDistributionWidget()
        self.basicWidget.changed.connect(self._onUIChanged)
        self.galaxyWidget.changed.connect(self._onUIChanged)
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

    def _onUIChanged(self):
        self._updateButtonStates()

    def _updateButtonStates(self):
        hasChanges = self.uiParameters != self.activeParameters
        self.reinitializationButton.setEnabled(hasChanges)

    def _reinitialize(self):
        self._onParamChanged()
        self.launchSimulation.emit(self.parameters)

    def _onDistributionTypeChanged(self, text):
        idx = 0 if text == "Basic" else 1
        self.distributionTab.setCurrentIndex(idx)
        self._onParamChanged()

    def _onParamChanged(self):
        pass

    def getParameters(self) -> SimulatorParameters:
        return self.parameters


class BasicDistributionWidget(QWidget):
    changed = pyqtSignal()

    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or BasicDistributionParameters()
        # PARAMETERS USER INTERFACE
        self.numParticles = QSpinBox()
        self.numParticles.setRange(100, 100000)
        self.numParticles.setValue(self.parameters.numParticles)
        self.numParticles.valueChanged.connect(self.changed.emit)
        self.positionScale = QDoubleSpinBox()
        self.positionScale.setRange(1.0, 100.0)
        self.positionScale.setValue(self.parameters.positionScale)
        self.positionScale.valueChanged.connect(self.changed.emit)
        self.velocityScale = QDoubleSpinBox()
        self.velocityScale.setRange(0.1, 10.0)
        self.velocityScale.setValue(self.parameters.velocityScale)
        self.velocityScale.valueChanged.connect(self.changed.emit)
        self.is3D = QCheckBox("3D Simulation")
        self.is3D.setChecked(self.parameters.is3D)
        self.is3D.stateChanged.connect(self.changed.emit)
        # MAIN LAYOUT
        group = QGroupBox("Basic Random Distribution Parameters")
        form = QFormLayout(group)
        form.addRow("Number of Particles:", self.numParticles)
        form.addRow("Position Scale:", self.positionScale)
        form.addRow("Velocity Scale:", self.velocityScale)
        form.addRow("", self.is3D)
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(group)

    def getParams(self):
        self.parameters.numParticles = self.numParticles.value()
        self.parameters.positionScale = self.positionScale.value()
        self.parameters.velocityScale = self.velocityScale.value()
        self.parameters.is3D = self.is3D.isChecked()
        return self.parameters



class GalaxyDistributionWidget(QWidget):
    changed = pyqtSignal()

    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or GalaxyDistributionParameters()
        # PARAMETERS USER INTERFACE
        self.numParticles = QSpinBox()
        self.numParticles.setRange(500, 100000)
        self.numParticles.setValue(self.parameters.numParticles)
        self.numParticles.valueChanged.connect(self.changed.emit)
        self.totalMass = QDoubleSpinBox()
        self.totalMass.setRange(100.0, 10000.0)
        self.totalMass.setValue(self.parameters.totalMass)
        self.totalMass.valueChanged.connect(self.changed.emit)
        self.radius = QDoubleSpinBox()
        self.radius.setRange(5.0, 100.0)
        self.radius.setValue(self.parameters.radius)
        self.radius.valueChanged.connect(self.changed.emit)
        self.height = QDoubleSpinBox()
        self.height.setRange(0.5, 30.0)
        self.height.setValue(self.parameters.height)
        self.height.valueChanged.connect(self.changed.emit)
        self.is3D = QCheckBox("3D Simulation")
        self.is3D.setChecked(self.parameters.is3D)
        self.is3D.stateChanged.connect(self.changed.emit)
        # MAIN LAYOUT
        group = QGroupBox("Galaxy Disk Distribution Parameters")
        form = QFormLayout(group)
        form.addRow("Number of Particles:", self.numParticles)
        form.addRow("Total Mass:", self.totalMass)
        form.addRow("Radius:", self.radius)
        form.addRow("Disk Height:", self.height)
        form.addRow("", self.is3D)
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(group)

    def getParams(self):
        self.parameters.numParticles = self.numParticles.value()
        self.parameters.totalMass = self.totalMass.value()
        self.parameters.radius = self.radius.value()
        self.parameters.height = self.height.value()
        self.parameters.is3D = self.is3D.isChecked()
        return self.parameters
