from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import *

from src.core.engine.parameters import SimulatorParameters, BasicDistributionParameters, GalaxyDistributionParameters


class SimulationConfigEditorDock(QDockWidget):
    launchSimulationPressed = pyqtSignal(object)
    resetSimulationPressed = pyqtSignal(object)

    def __init__(self, initialParameters: SimulatorParameters = None, parent=None):
        super().__init__("Simulation Configuration", parent)
        self.activeParameters = initialParameters or SimulatorParameters()
        self.uiParameters = self.activeParameters
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
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
        self.distributionEditor = QStackedWidget()
        self.basicWidget = BasicDistributionWidget()
        self.galaxyWidget = GalaxyDistributionWidget()
        self.basicWidget.changed.connect(self._onUIChanged)
        self.galaxyWidget.changed.connect(self._onUIChanged)
        self.distributionEditor.addWidget(self.basicWidget)
        self.distributionEditor.addWidget(self.galaxyWidget)

        # MAIN LAYOUT & CONTAINER
        self.container = QWidget()
        self.launchButton = QPushButton("Launch Simulation")
        self.launchButton.clicked.connect(self._launchSimulation)
        self.resetButton = QPushButton("Reset Simulation")
        self.resetButton.clicked.connect(self._resetSimulation)
        self.resetButton.setEnabled(False)
        buttonsLayout = QHBoxLayout()
        buttonsLayout.addWidget(self.launchButton)
        buttonsLayout.addWidget(self.resetButton)
        mainLayout = QVBoxLayout(self.container)
        mainLayout.addWidget(self.simulationTab)
        mainLayout.addWidget(self.distributionEditor)
        mainLayout.addStretch()
        mainLayout.addLayout(buttonsLayout)
        self.setWidget(self.container)
        self.hasRunBefore = False
        self._initializeUi()
        self._updateButtonStates()

    def _initializeUi(self):
        p = self.activeParameters
        self.nameEdit.blockSignals(True)
        self.distributionTypeCombo.blockSignals(True)
        self.calculatorTypeCombo.blockSignals(True)
        self.integratorTypeCombo.blockSignals(True)
        self.gravitationalConstantSpin.blockSignals(True)
        self.timeStepSpin.blockSignals(True)
        self.thetaValueSpin.blockSignals(True)
        self.endlessCheck.blockSignals(True)
        self.maxTimeSpin.blockSignals(True)
        self.saveResultsCheck.blockSignals(True)
        try:
            self.nameEdit.setText(p.name)
            self.distributionTypeCombo.setCurrentText(p.distributionType.title())
            self.calculatorTypeCombo.setCurrentText(p.calculatorType.replace("_", "-"))
            self.integratorTypeCombo.setCurrentText(p.integratorType)
            self.gravitationalConstantSpin.setValue(p.gravitationalConstant)
            self.timeStepSpin.setValue(p.timeStep)
            self.thetaValueSpin.setValue(p.theta)
            self.endlessCheck.setChecked(p.endless)
            self.maxTimeSpin.setValue(p.maxTime)
            self.saveResultsCheck.setChecked(p.saveResults)
            idx = 0 if p.distributionType.upper() == "BASIC" else 1
            self.distributionEditor.setCurrentIndex(idx)
            self.galaxyWidget.syncUIFromParameters(p.galaxyDistributionParameters)
            self.basicWidget.syncUIFromParameters(p.basicDistributionParameters)
            self.uiParameters = p
        finally:
            # Re-enable signals
            self.nameEdit.blockSignals(False)
            self.distributionTypeCombo.blockSignals(False)
            self.calculatorTypeCombo.blockSignals(False)
            self.integratorTypeCombo.blockSignals(False)
            self.gravitationalConstantSpin.blockSignals(False)
            self.timeStepSpin.blockSignals(False)
            self.thetaValueSpin.blockSignals(False)
            self.endlessCheck.blockSignals(False)
            self.maxTimeSpin.blockSignals(False)
            self.saveResultsCheck.blockSignals(False)

    def _onUIChanged(self):
        self.uiParameters.name = self.nameEdit.text()
        self.uiParameters.distributionType = self.distributionTypeCombo.currentText().upper()
        self.uiParameters.calculatorType = self.calculatorTypeCombo.currentText().replace("-", "_").upper()
        self.uiParameters.integratorType = self.integratorTypeCombo.currentText()
        self.uiParameters.gravitationalConstant = self.gravitationalConstantSpin.value()
        self.uiParameters.timeStep = self.timeStepSpin.value()
        self.uiParameters.theta = self.thetaValueSpin.value()
        self.uiParameters.endless = self.endlessCheck.isChecked()
        self.uiParameters.maxTime = self.maxTimeSpin.value()
        self.uiParameters.saveResults = self.saveResultsCheck.isChecked()
        if self.uiParameters.distributionType == "GALAXY":
            self.uiParameters.galaxyDistributionParameters = self.galaxyWidget.getParameters()
        else:
            self.uiParameters.basicDistributionParameters = self.basicWidget.getParameters()
        self._updateButtonStates()

    def _updateButtonStates(self):
        hasChanges = self.getUiParameters() != self.activeParameters
        self.launchButton.setEnabled(hasChanges or not self.hasRunBefore)
        self.resetButton.setEnabled(self.hasRunBefore)

    def _launchSimulation(self):
        self._onUIChanged()
        self.activeParameters = self.uiParameters
        self.hasRunBefore = True
        self.launchSimulationPressed.emit(self.activeParameters)
        self._updateButtonStates()

    def _resetSimulation(self):
        if self.hasRunBefore:
            self.resetSimulationPressed.emit(self.activeParameters)

    def _onDistributionTypeChanged(self, text):
        idx = 0 if text == "Basic" else 1
        self.distributionEditor.setCurrentIndex(idx)
        self._onUIChanged()

    def getUiParameters(self):
        return self.uiParameters

    def getParameters(self):
        return self.activeParameters

    def notifySimulationStarted(self):
        self.hasRunBefore = True
        self._updateButtonStates()


class BasicDistributionWidget(QWidget):
    changed = pyqtSignal()

    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or BasicDistributionParameters()
        # PARAMETERS USER INTERFACE
        self.numParticles = QSpinBox()
        self.numParticles.setRange(100, 100000)
        self.numParticles.setValue(self.parameters.nbParticles)
        self.numParticles.valueChanged.connect(self.changed.emit)
        self.positionScale = QDoubleSpinBox()
        self.positionScale.setRange(1.0, 100.0)
        self.positionScale.setValue(self.parameters.positionScale)
        self.positionScale.valueChanged.connect(self.changed.emit)
        self.velocityScale = QDoubleSpinBox()
        self.velocityScale.setRange(0.1, 10.0)
        self.velocityScale.setValue(self.parameters.velocityScale)
        self.velocityScale.valueChanged.connect(self.changed.emit)
        # MAIN LAYOUT
        group = QGroupBox("Basic Random Distribution Parameters")
        form = QFormLayout(group)
        form.addRow("Number of Particles:", self.numParticles)
        form.addRow("Position Scale:", self.positionScale)
        form.addRow("Velocity Scale:", self.velocityScale)
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(group)

    def getParameters(self):
        self.parameters.nbParticles = self.numParticles.value()
        self.parameters.positionScale = self.positionScale.value()
        self.parameters.velocityScale = self.velocityScale.value()
        return self.parameters

    def syncUIFromParameters(self, parameters: BasicDistributionParameters):
        self.numParticles.blockSignals(True)
        self.positionScale.blockSignals(True)
        self.velocityScale.blockSignals(True)
        try:
            self.numParticles.setValue(parameters.nbParticles)
            self.positionScale.setValue(parameters.positionScale)
            self.velocityScale.setValue(parameters.velocityScale)
            self.parameters = parameters
        finally:
            self.numParticles.blockSignals(False)
            self.positionScale.blockSignals(False)
            self.velocityScale.blockSignals(False)



class GalaxyDistributionWidget(QWidget):
    changed = pyqtSignal()

    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or GalaxyDistributionParameters()
        # PARAMETERS USER INTERFACE
        self.numParticles = QSpinBox()
        self.numParticles.setRange(500, 100000)
        self.numParticles.setValue(self.parameters.nbParticles)
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
        # MAIN LAYOUT
        group = QGroupBox("Galaxy Disk Distribution Parameters")
        form = QFormLayout(group)
        form.addRow("Number of Particles:", self.numParticles)
        form.addRow("Total Mass:", self.totalMass)
        form.addRow("Radius:", self.radius)
        form.addRow("Disk Height:", self.height)
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(group)

    def getParameters(self):
        self.parameters.nbParticles = self.numParticles.value()
        self.parameters.totalMass = self.totalMass.value()
        self.parameters.radius = self.radius.value()
        self.parameters.height = self.height.value()
        return self.parameters

    def syncUIFromParameters(self, parameters: GalaxyDistributionParameters):
        self.numParticles.blockSignals(True)
        self.totalMass.blockSignals(True)
        self.radius.blockSignals(True)
        self.height.blockSignals(True)
        try:
            self.numParticles.setValue(parameters.nbParticles)
            self.totalMass.setValue(parameters.totalMass)
            self.radius.setValue(parameters.radius)
            self.height.setValue(parameters.height)
            self.parameters = parameters
        finally:
            self.numParticles.blockSignals(False)
            self.totalMass.blockSignals(False)
            self.radius.blockSignals(False)
            self.height.blockSignals(False)
            self.is3D.blockSignals(False)
