import copy
from random import randint
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import *

from src.gui.parameters import SimulatorParameters, BasicDistributionParameters, GalaxyDistributionParameters


class SimulationConfigEditorDock(QDockWidget):
    launchSimulationPressed = pyqtSignal(object)
    resetSimulationPressed = pyqtSignal(object)

    def __init__(self, initialParameters: SimulatorParameters, parent=None):
        super().__init__("Simulation Configuration", parent)
        self.activeParameters = copy.deepcopy(initialParameters or SimulatorParameters())
        self.uiParameters = copy.deepcopy(self.activeParameters)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        self.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetClosable)
        self.setFloating(False)
        self.hasRunBefore = False

        # MAIN SIMULATION PARAMETERS
        self.simulationTab = QWidget()
        self.nameEdit = QLineEdit("My N-Body Simulation")
        self.nameEdit.textChanged.connect(self._onUIChanged)
        self.saveResultsCheck = QCheckBox("Save Results to Disk")
        self.saveResultsCheck.stateChanged.connect(self._onUIChanged)
        simulationInfoGroup = QGroupBox("Simulation Information")
        simulationInfoForm = QFormLayout(simulationInfoGroup)
        simulationInfoForm.addRow("Name:", self.nameEdit)
        simulationInfoForm.addRow("", self.saveResultsCheck)
        self.distributionTypeCombo = QComboBox()
        self.distributionTypeCombo.addItem('Galaxy', userData='GALAXY')
        self.distributionTypeCombo.addItem('Basic Random', userData='BASIC')
        self.distributionTypeCombo.currentIndexChanged.connect(self._onDistributionTypeChanged)
        self.seedEdit = QLineEdit()
        self.seedEdit.setValidator(QIntValidator(0, 2 ** 31 - 1))
        self.seedEdit.editingFinished.connect(self._onUIChanged)
        self.randomSeedButton = QPushButton("Randomize")
        self.randomSeedButton.clicked.connect(self._randomizeSeed)
        self.calculatorTypeCombo = QComboBox()
        self.calculatorTypeCombo.addItem('Newton', userData='NEWTON')
        self.calculatorTypeCombo.addItem('Barnes-Hut', userData='BARNES_HUT')
        self.calculatorTypeCombo.currentIndexChanged.connect(self._onUIChanged)
        self.integratorTypeCombo = QComboBox()
        self.integratorTypeCombo.addItem('RK4', userData='RK4')
        self.integratorTypeCombo.addItem('Euler Explicit', userData='EULER_EXPLICIT')
        self.integratorTypeCombo.currentIndexChanged.connect(self._onUIChanged)
        self.deviceCombo = QComboBox()
        self.deviceCombo.addItem('GPU (CUDA)', userData='GPU')
        self.deviceCombo.addItem('CPU', userData='CPU')
        self.deviceCombo.currentIndexChanged.connect(self._onUIChanged)
        self.numParticles = QSpinBox()
        self.numParticles.setRange(500, 1000000)
        self.numParticles.valueChanged.connect(self._onUIChanged)
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
        seedLayout = QHBoxLayout()
        seedLayout.addWidget(self.seedEdit)
        seedLayout.addWidget(self.randomSeedButton)
        simulationCoreGroup = QGroupBox("Core Simulation Parameters")
        simulationCoreForm = QFormLayout(simulationCoreGroup)
        simulationCoreForm.addRow("Distribution Type:", self.distributionTypeCombo)
        simulationCoreForm.addRow("Seed:", seedLayout)
        simulationCoreForm.addRow("Force Calculator:", self.calculatorTypeCombo)
        simulationCoreForm.addRow("Integrator:", self.integratorTypeCombo)
        simulationCoreForm.addRow("Device:", self.deviceCombo)
        simulationCoreForm.addRow("Number of Particles:", self.numParticles)
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
        parameters = copy.deepcopy(self.activeParameters)
        self.nameEdit.blockSignals(True)
        self.distributionTypeCombo.blockSignals(True)
        self.seedEdit.blockSignals(True)
        self.calculatorTypeCombo.blockSignals(True)
        self.integratorTypeCombo.blockSignals(True)
        self.deviceCombo.blockSignals(True)
        self.gravitationalConstantSpin.blockSignals(True)
        self.timeStepSpin.blockSignals(True)
        self.thetaValueSpin.blockSignals(True)
        self.endlessCheck.blockSignals(True)
        self.maxTimeSpin.blockSignals(True)
        self.saveResultsCheck.blockSignals(True)
        try:
            self.nameEdit.setText(parameters.name)
            self.seedEdit.setText(str(randint(0, 2 ** 31 - 1)) if parameters.seed is None else str(parameters.seed))
            self._setComboData(self.distributionTypeCombo, parameters.distributionType)
            self._setComboData(self.calculatorTypeCombo, parameters.calculatorType)
            self._setComboData(self.integratorTypeCombo, parameters.integratorType)
            self._setComboData(self.deviceCombo, parameters.device)
            self.gravitationalConstantSpin.setValue(parameters.gravitationalConstant)
            self.timeStepSpin.setValue(parameters.timeStep)
            self.thetaValueSpin.setValue(parameters.theta)
            self.endlessCheck.setChecked(parameters.endless)
            self.maxTimeSpin.setValue(parameters.maxTime)
            self.saveResultsCheck.setChecked(parameters.saveResults)
            idx = 0 if parameters.distributionType.upper() == "BASIC" else 1
            self.distributionEditor.setCurrentIndex(idx)
            self.galaxyWidget.syncUIFromParameters(parameters.galaxyDistributionParameters)
            self.basicWidget.syncUIFromParameters(parameters.basicDistributionParameters)
            self.uiParameters = copy.deepcopy(parameters)
        finally:
            self.nameEdit.blockSignals(False)
            self.distributionTypeCombo.blockSignals(False)
            self.seedEdit.blockSignals(False)
            self.calculatorTypeCombo.blockSignals(False)
            self.integratorTypeCombo.blockSignals(False)
            self.deviceCombo.blockSignals(False)
            self.gravitationalConstantSpin.blockSignals(False)
            self.timeStepSpin.blockSignals(False)
            self.thetaValueSpin.blockSignals(False)
            self.endlessCheck.blockSignals(False)
            self.maxTimeSpin.blockSignals(False)
            self.saveResultsCheck.blockSignals(False)

    def _onUIChanged(self):
        self.uiParameters.name = self.nameEdit.text()
        self.uiParameters.distributionType = self.distributionTypeCombo.currentData()
        seedText = self.seedEdit.text().strip()
        self.uiParameters.seed = int(seedText) if seedText else None
        self.uiParameters.calculatorType = self.calculatorTypeCombo.currentData()
        self.uiParameters.integratorType = self.integratorTypeCombo.currentData()
        self.uiParameters.device = self.deviceCombo.currentData()
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
        self.activeParameters = copy.deepcopy(self.uiParameters)
        self.hasRunBefore = True
        self.launchSimulationPressed.emit(copy.deepcopy(self.activeParameters))
        self._updateButtonStates()

    def _resetSimulation(self):
        if self.hasRunBefore:
            self.resetSimulationPressed.emit(copy.deepcopy(self.activeParameters))

    def _randomizeSeed(self):
        seed = randint(0, 2 ** 31 - 1)
        self.seedEdit.setText(str(seed))
        self._onUIChanged()

    def _onDistributionTypeChanged(self, text):
        idx = 0 if text == "Basic" else 1
        self.distributionEditor.setCurrentIndex(idx)
        self._onUIChanged()

    def getUiParameters(self):
        return copy.deepcopy(self.uiParameters)

    def getParameters(self):
        return copy.deepcopy(self.activeParameters)

    @staticmethod
    def _setComboData(combo: QComboBox, value):
        index = combo.findData(value)
        if index != -1:
            combo.setCurrentIndex(index)


class BasicDistributionWidget(QWidget):
    changed = pyqtSignal()

    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or BasicDistributionParameters()
        # PARAMETERS USER INTERFACE
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
        form.addRow("Position Scale:", self.positionScale)
        form.addRow("Velocity Scale:", self.velocityScale)
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(group)

    def getParameters(self):
        self.parameters.positionScale = self.positionScale.value()
        self.parameters.velocityScale = self.velocityScale.value()
        return self.parameters

    def syncUIFromParameters(self, parameters: BasicDistributionParameters):
        self.positionScale.blockSignals(True)
        self.velocityScale.blockSignals(True)
        try:
            self.positionScale.setValue(parameters.positionScale)
            self.velocityScale.setValue(parameters.velocityScale)
            self.parameters = parameters
        finally:
            self.positionScale.blockSignals(False)
            self.velocityScale.blockSignals(False)



class GalaxyDistributionWidget(QWidget):
    changed = pyqtSignal()

    def __init__(self, parameters=None):
        super().__init__()
        self.parameters = parameters or GalaxyDistributionParameters()
        # PARAMETERS USER INTERFACE
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
        self.bulgeFraction = QDoubleSpinBox()
        self.bulgeFraction.setRange(0.0, 1.0)
        self.bulgeFraction.setSingleStep(0.05)
        self.bulgeFraction.setDecimals(2)
        self.bulgeFraction.setValue(self.parameters.bulgeFraction)
        self.bulgeFraction.valueChanged.connect(self._updateFractionLimits)
        self.haloFraction = QDoubleSpinBox()
        self.haloFraction.setRange(0.0, 1.0)
        self.haloFraction.setSingleStep(0.05)
        self.haloFraction.setDecimals(2)
        self.haloFraction.setValue(self.parameters.haloFraction)
        self.haloFraction.valueChanged.connect(self._updateFractionLimits)
        self.plummerRadius = QDoubleSpinBox()
        self.plummerRadius.setRange(0.1, 1000.0)
        self.plummerRadius.setDecimals(2)
        self.plummerRadius.setValue(self.parameters.plummerRadius)
        self.plummerRadius.valueChanged.connect(self.changed.emit)
        self.haloRadius = QDoubleSpinBox()
        self.haloRadius.setRange(1.0, 1000.0)
        self.haloRadius.setDecimals(2)
        self.haloRadius.setValue(self.parameters.haloRadius)
        self.haloRadius.valueChanged.connect(self.changed.emit)
        self.velocityDispersion = QDoubleSpinBox()
        self.velocityDispersion.setRange(0.0, 10.0)
        self.velocityDispersion.setDecimals(3)
        self.velocityDispersion.setSingleStep(0.01)
        self.velocityDispersion.setValue(self.parameters.velocityDispersion)
        self.velocityDispersion.valueChanged.connect(self.changed.emit)
        # MAIN LAYOUT
        group = QGroupBox("Galaxy Distribution Parameters")
        diskGroup = QGroupBox("Disk")
        diskForm = QFormLayout(diskGroup)
        diskForm.addRow("Total Mass:", self.totalMass)
        diskForm.addRow("Radius:", self.radius)
        diskForm.addRow("Disk Height:", self.height)
        structureGroup = QGroupBox("Galaxy Structure")
        structureForm = QFormLayout(structureGroup)
        structureForm.addRow("Bulge Fraction:", self.bulgeFraction)
        structureForm.addRow("Halo Fraction:", self.haloFraction)
        structureForm.addRow("Plummer Radius:", self.plummerRadius)
        structureForm.addRow("Halo Radius:", self.haloRadius)
        velocityGroup = QGroupBox("Velocity")
        velocityForm = QFormLayout(velocityGroup)
        velocityForm.addRow("Velocity Dispersion:", self.velocityDispersion)
        groupLayout = QVBoxLayout(group)
        groupLayout.addWidget(diskGroup)
        groupLayout.addWidget(structureGroup)
        groupLayout.addWidget(velocityGroup)
        mainLayout = QVBoxLayout(self)
        mainLayout.addWidget(group)

    def getParameters(self):
        self.parameters.totalMass = self.totalMass.value()
        self.parameters.radius = self.radius.value()
        self.parameters.height = self.height.value()
        self.parameters.bulgeFraction = self.bulgeFraction.value()
        self.parameters.haloFraction = self.haloFraction.value()
        self.parameters.velocityDispersion = self.velocityDispersion.value()
        self.parameters.plummerRadius = self.plummerRadius.value()
        self.parameters.haloRadius = self.haloRadius.value()
        return self.parameters

    def syncUIFromParameters(self, parameters: GalaxyDistributionParameters):
        self.totalMass.blockSignals(True)
        self.radius.blockSignals(True)
        self.height.blockSignals(True)
        self.bulgeFraction.blockSignals(True)
        self.haloFraction.blockSignals(True)
        self.velocityDispersion.blockSignals(True)
        self.plummerRadius.blockSignals(True)
        self.haloRadius.blockSignals(True)
        try:
            self.totalMass.setValue(parameters.totalMass)
            self.radius.setValue(parameters.radius)
            self.height.setValue(parameters.height)
            self.bulgeFraction.setValue(parameters.bulgeFraction)
            self.haloFraction.setValue(parameters.haloFraction)
            self.velocityDispersion.setValue(parameters.velocityDispersion)
            self.plummerRadius.setValue(parameters.plummerRadius)
            self.haloRadius.setValue(parameters.haloRadius)
            self.parameters = parameters
        finally:
            self.totalMass.blockSignals(False)
            self.radius.blockSignals(False)
            self.height.blockSignals(False)
            self.bulgeFraction.blockSignals(False)
            self.haloFraction.blockSignals(False)
            self.velocityDispersion.blockSignals(False)
            self.plummerRadius.blockSignals(False)
            self.haloRadius.blockSignals(False)

    def _updateFractionLimits(self):
        bulgeFraction = self.bulgeFraction.value()
        haloFraction = self.haloFraction.value()
        self.bulgeFraction.blockSignals(True)
        self.haloFraction.blockSignals(True)
        self.bulgeFraction.setMaximum(1.0 - haloFraction)
        self.haloFraction.setMaximum(1.0 - bulgeFraction)
        self.bulgeFraction.blockSignals(False)
        self.haloFraction.blockSignals(False)
        self.changed.emit()
