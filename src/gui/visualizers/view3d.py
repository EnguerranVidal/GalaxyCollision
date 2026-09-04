from dataclasses import dataclass

import numpy as np
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GL import *

from PyQt5.QtCore import Qt, pyqtSignal, QTimer
from PyQt5.QtWidgets import *

from src.gui.visualizers.renderers import ParticlesRenderer, SinglePointRenderer, VectorFieldRenderer, GridRenderer
from src.gui.settings import ViewSettings
from src.gui.solver.general import State


class Camera:
    def __init__(self, minimumZoom: float = 1.25, maximumZoom: float = 5000.0):
        self.zoom = 5.0
        self.rotationX = 45.0
        self.rotationY = 225.0
        self.minimumZoom = float(minimumZoom)
        self.maximumZoom = float(maximumZoom)

    def apply(self):
        glLoadIdentity()
        glTranslatef(0, 0, -self.zoom)
        glRotatef(self.rotationX, 1, 0, 0)
        glRotatef(self.rotationY, 0, 1, 0)
        glRotatef(-90, 1, 0, 0)

    def rotate(self, dx, dy):
        self.rotationY += dx * 0.5
        self.rotationX += dy * 0.5
        self.rotationX = max(-89.0, min(89.0, self.rotationX))

    def zoomBy(self, delta):
        self.zoom *= (0.9 if delta > 0 else 1.1)
        self.zoom = max(self.minimumZoom, min(self.maximumZoom, self.zoom))

    def getPosition(self):
        rotationX = np.deg2rad(self.rotationX)
        rotationY = np.deg2rad(self.rotationY)
        rotationFix = np.deg2rad(-90)
        camera = np.array([0, 0, self.zoom], dtype=float)
        xRotation = np.array([[1, 0, 0], [0, np.cos(-rotationX), -np.sin(-rotationX)], [0, np.sin(-rotationX), np.cos(-rotationX)]])
        yRotation = np.array([[np.cos(-rotationY), 0, np.sin(-rotationY)], [0, 1, 0], [-np.sin(-rotationY), 0, np.cos(-rotationY)]])
        fixRotation = np.array([[1, 0, 0], [0, np.cos(-rotationFix), -np.sin(-rotationFix)], [0, np.sin(-rotationFix), np.cos(-rotationFix)]])
        return fixRotation @ (yRotation @ (xRotation @ camera))

    def setZoomLimits(self, minimumZoom: float, maximumZoom: float):
        self.minimumZoom = float(minimumZoom)
        self.maximumZoom = float(maximumZoom)
        self.zoom = max(self.minimumZoom, min(self.maximumZoom, self.zoom))


@dataclass
class SelectedParticle:
    particleIndex: int | None = None
    particleGroup: str | None = None

    def unselect(self):
        self.particleIndex = None
        self.particleGroup = None

    def select(self, particleIndex, particleGroup):
        self.particleIndex = particleIndex
        self.particleGroup = particleGroup


class Universe3dViewWidget(QOpenGLWidget):
    cameraChanged = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        glutInit()
        self.setMouseTracking(True)
        self.lastPosX, self.lastPosY = 0, 0
        self.massCenter = np.zeros(3, dtype=np.float32)
        self.viewSettings = ViewSettings()
        self.camera = Camera(self.viewSettings.minimumExtent, self.viewSettings.maximumExtent)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._showContextMenu)

        # PARTICLE SELECTION
        self.selectedParticle = SelectedParticle()
        self.lastPositions = {}
        self._pressX, self._pressY = 0, 0
        self._dragged = False
        self._ignoreNextRelease = False
        self.focusOnSelectedParticle = False
        self.pulsePhase = 0.0
        self.pulseTimer = QTimer(self)
        self.pulseTimer.setInterval(16)
        self.pulseTimer.timeout.connect(self._onPulseTick)

        # RENDERERS AND OBJECT BUFFERS
        self.groupColors = {}
        self.pendingObjectBufferUpdates = {}
        self.gridRenderer = GridRenderer(self.viewSettings.minimumExtent, self.viewSettings.maximumExtent)
        self.particlesRenderer = ParticlesRenderer()
        self.barycenterRenderer = SinglePointRenderer()
        self.selectionRenderer = SinglePointRenderer()
        self.velocityVectorRenderer = VectorFieldRenderer(color=(0.35, 0.9, 1.0, 0.9))
        self.accelerationVectorRenderer = VectorFieldRenderer(color=(1.0, 0.45, 0.2, 0.9))
        self.velocityVectorVertices = None
        self.accelerationVectorVertices = None

        # TIME VALUE OVERLAY
        self.timeOverlay = QLabel(self)
        self.timeOverlay.setText("Time: 0.000")
        self.timeOverlay.setAttribute(Qt.WA_TransparentForMouseEvents, True)
        self.timeOverlay.adjustSize()
        self._placeTimeOverlay()

    def initializeGL(self):
        glClearColor(0, 0, 0, 1.0)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LEQUAL)
        glClearDepth(1.0)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_MULTISAMPLE)
        glEnable(GL_LINE_SMOOTH)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
        glEnable(GL_POINT_SMOOTH)
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST)
        glEnable(GL_POINT_SPRITE)
        glEnable(GL_PROGRAM_POINT_SIZE)
        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
        glShadeModel(GL_SMOOTH)
        self.particlesRenderer.initialize()
        self.barycenterRenderer.initialize()
        self.selectionRenderer.initialize()
        self.velocityVectorRenderer.initialize()
        self.accelerationVectorRenderer.initialize()
        self._updateProjection(max(self.width(), 1), max(self.height(), 1))

    def paintGL(self):
        self._uploadPendingObjectBuffers()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.camera.apply()
        plotOrigin = self._plotOrigin()
        # FOCUS ON SELECTED PARTICLE
        if self.focusOnSelectedParticle:
            particlePosition = self._selectedWorldPosition()
            if particlePosition is not None:
                glTranslatef(-float(particlePosition[0]- plotOrigin[0]), -float(particlePosition[1] - plotOrigin[1]), -float(particlePosition[2] - plotOrigin[2]))
        self.gridRenderer.render(self.camera.zoom)
        glTranslatef(-float(plotOrigin[0]),  -float(plotOrigin[1]), -float(plotOrigin[2]))
        glDisable(GL_LIGHTING)
        if self.viewSettings.showBarycenter:
            self.barycenterRenderer.render(self.massCenter, pointSize=14.0)
        self.particlesRenderer.renderAll(pointSize=4.0, refDistance=max(self.camera.zoom, 0.5), minSize=1.5, maxSize=16.0)
        if self.viewSettings.showVelocityVectors and self.velocityVectorVertices is not None:
            self.velocityVectorRenderer.update(self.velocityVectorVertices)
            self.velocityVectorRenderer.render(lineWidth=1.5)
        if self.viewSettings.showAccelerationVectors and self.accelerationVectorVertices is not None:
            self.accelerationVectorRenderer.update(self.accelerationVectorVertices)
            self.accelerationVectorRenderer.render(lineWidth=1.5)
        if self.selectedParticle.particleIndex is not None:
            particlePosition = self._selectedWorldPosition()
            if particlePosition is not None:
                drawPos = np.asarray(particlePosition, dtype=np.float32)
                pulse = 1.0 + 0.35 * np.sin(self.pulsePhase)
                self.selectionRenderer.render(drawPos, pointSize=20.0 * pulse, color=(1.0, 0.95, 0.2, 1.0))

    def updateState(self, state: State):
        self.pendingObjectBufferUpdates.clear()
        # PARTICLES POSITIONS
        for groupIndex, positionArray in state.positions.items():
            if groupIndex not in self.groupColors:
                self.groupColors[groupIndex] = (0.9, 0.7, 0.3, 0.95)
            self.pendingObjectBufferUpdates[groupIndex] = positionArray
        self.massCenter = np.asarray(state.massCenter, dtype=np.float32).reshape(3)
        # SELECTION UPDATE
        self.lastPositions = {groupIndex: np.ascontiguousarray(positionArray, dtype=np.float32) for groupIndex, positionArray in state.positions.items()}
        if self.selectedParticle.particleIndex is not None:
            groupName = self.selectedParticle.particleGroup
            particleIndex = self.selectedParticle.particleIndex
            positions = self.lastPositions.get("default" if groupName is None else groupName, None)
            if positions is None or particleIndex >= len(positions):
                self.selectedParticle.unselect()
                self.focusOnSelectedParticle = False
                self.pulseTimer.stop()
        # VELOCITY ARRAY FIELD
        if self.viewSettings.showVelocityVectors and state.velocities:
            chunks = []
            nbTotal = sum(len(pos) for pos in state.positions.values())
            subSample = 1 if nbTotal <= 8000 else max(1, nbTotal // 8000)
            for groupIndex, positions in state.positions.items():
                velocities = state.velocities.get(groupIndex)
                if velocities is None or len(velocities) != len(positions):
                    continue
                vectorLength, reference = self.viewSettings.velocityVectorLength, self.viewSettings.referenceVelocity
                vectors = VectorFieldRenderer.buildVectors(positions, velocities, vectorLength, reference, subsample=subSample)
                chunks.append(vectors)
            self.velocityVectorVertices = np.vstack(chunks) if chunks else None
        else:
            self.velocityVectorVertices = None
        # ACCELERATION ARRAY FIELD
        if self.viewSettings.showAccelerationVectors and state.accelerations:
            chunks = []
            nbTotal = sum(len(pos) for pos in state.positions.values())
            subsample = 1 if nbTotal <= 8000 else max(1, nbTotal // 8000)
            for groupIndex, positions in state.positions.items():
                accelerations = state.accelerations.get(groupIndex)
                if accelerations is None or len(accelerations) != len(positions):
                    continue
                vectorLength, reference = self.viewSettings.accelerationVectorLength, self.viewSettings.referenceAcceleration
                vectors = VectorFieldRenderer.buildVectors(positions, accelerations, vectorLength, reference, subsample=subsample)
                chunks.append(vectors)
            self.accelerationVectorVertices = np.vstack(chunks) if chunks else None
        else:
            self.accelerationVectorVertices = None
        # SIMULATION TIME OVERLAY
        self.timeOverlay.setText(f"Time: {state.time:.3f}")
        self.timeOverlay.adjustSize()
        self._placeTimeOverlay()
        self.update()

    def _applyExtentLimits(self):
        self.camera.setZoomLimits(self.viewSettings.minimumExtent, self.viewSettings.maximumExtent)
        self.gridRenderer.setExtentLimits(self.viewSettings.minimumExtent, self.viewSettings.maximumExtent)
        if self.isValid():
            self.makeCurrent()
            self._updateProjection(max(self.width(), 1), max(self.height(), 1))
            self.doneCurrent()
        self.update()

    def setViewSettings(self, viewSettings: ViewSettings):
        self.viewSettings = viewSettings
        if self.viewSettings.maximumExtent < self.viewSettings.minimumExtent:
            self.viewSettings.maximumExtent = self.viewSettings.minimumExtent
        self._applyExtentLimits()
        if not self.viewSettings.showBarycenter and "barycenter" in self.particlesRenderer.counts:
            self.particlesRenderer.counts["barycenter"] = 0
        self.update()

    def _placeTimeOverlay(self):
        margin = 12
        x, y = self.width() - self.timeOverlay.width() - margin, self.height() - self.timeOverlay.height() - margin
        self.timeOverlay.move(max(0, x), max(0, y))
        self.timeOverlay.raise_()

    def _onPulseTick(self):
        if self.selectedParticle.particleIndex is None:
            self.pulseTimer.stop()
            return
        self.pulsePhase += 0.15
        self.update()

    def _uploadPendingObjectBuffers(self):
        if not self.pendingObjectBufferUpdates:
            return
        for groupIndex, positions in self.pendingObjectBufferUpdates.items():
            if groupIndex not in self.particlesRenderer.colors:
                self.particlesRenderer.createGroup(groupIndex, self.groupColors.get(groupIndex, (1.0, 1.0, 1.0, 0.9)))
            self.particlesRenderer.updateGroupPositions(groupIndex, positions)
        self.pendingObjectBufferUpdates.clear()

    def _pickParticle(self, xMouse: int, yMouse: int, maxPixelDistance: float = 15.0):
        if not self.lastPositions:
            return None, None
        self.makeCurrent()
        ratio = float(self.devicePixelRatioF())
        xMouse = xMouse * ratio
        yMouse = yMouse * ratio
        maxPixelDistance = maxPixelDistance * ratio
        viewModel = (GLdouble * 16)()
        viewProjection = (GLdouble * 16)()
        viewPort = (GLint * 4)()
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        self.camera.apply()
        plotOrigin = self._plotOrigin()
        glTranslatef( -float(plotOrigin[0]), -float(plotOrigin[1]), -float(plotOrigin[2]))
        glGetDoublev(GL_MODELVIEW_MATRIX, viewModel)
        glGetDoublev(GL_PROJECTION_MATRIX, viewProjection)
        glGetIntegerv(GL_VIEWPORT, viewPort)
        glPopMatrix()
        bestDistanceSquared = maxPixelDistance * maxPixelDistance
        bestGroup, bestIndex = None, None
        yView = viewPort[3] - yMouse
        for groupIndex, positions in self.lastPositions.items():
            for i, (x, y, z) in enumerate(positions):
                xWindow, yWindow, zWindow = gluProject(float(x), float(y), float(z), viewModel, viewProjection, viewPort)
                if zWindow < 0.0 or zWindow > 1.0:
                    continue
                xDelta, yDelta = xWindow - xMouse, yWindow - yView
                squaredDistance = xDelta * xDelta + yDelta * yDelta
                if squaredDistance < bestDistanceSquared:
                    bestDistanceSquared = squaredDistance
                    bestGroup, bestIndex = groupIndex, i
        return bestGroup, bestIndex

    def _selectedWorldPosition(self):
        groupName = self.selectedParticle.particleGroup
        particleIndex = self.selectedParticle.particleIndex
        if groupName is None or particleIndex is None:
            return None
        positions = self.lastPositions.get(groupName)
        if positions is None or particleIndex >= len(positions):
            return None
        return positions[particleIndex]

    def _plotOrigin(self):
        return self.massCenter if self.viewSettings.centerOnBarycenter else np.zeros(3, dtype=np.float32)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._pressX = event.x()
            self._pressY = event.y()
            self.lastPosX = event.x()
            self.lastPosY = event.y()
            self._dragged = False

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.LeftButton:
            dx = event.x() - self.lastPosX
            dy = event.y() - self.lastPosY
            if abs(event.x() - self._pressX) > 5 or abs(event.y() - self._pressY) > 5:
                self._dragged = True
            self.camera.rotate(dx, dy)
            self.lastPosX = event.x()
            self.lastPosY = event.y()
            self.update()
            self.cameraChanged.emit()

    def mouseReleaseEvent(self, event):
        if event.button() != Qt.LeftButton or self._dragged:
            return
        if self._ignoreNextRelease:
            self._ignoreNextRelease = False
            return
        group, index = self._pickParticle(event.x(), event.y())
        if group is None:
            self.selectedParticle.unselect()
            self.focusOnSelectedParticle = False
            self.pulseTimer.stop()
        else:
            self.selectedParticle.select(index, group)
            self.focusOnSelectedParticle = False
            self.pulsePhase = 0.0
            self.pulseTimer.start()
        self.update()

    def mouseDoubleClickEvent(self, event):
        if event.button() != Qt.LeftButton:
            return
        group, index = self._pickParticle(event.x(), event.y())
        if group is None:
            return
        self.selectedParticle.select(index, group)
        self.focusOnSelectedParticle = True
        self._ignoreNextRelease = True
        self.pulsePhase = 0.0
        self.pulseTimer.start()
        self.update()

    def wheelEvent(self, event):
        delta = event.angleDelta().y() / 120.0
        self.camera.zoomBy(delta)
        self.update()
        self.cameraChanged.emit()

    def resetView(self):
        self.camera = Camera(self.viewSettings.minimumExtent, self.viewSettings.maximumExtent)
        self.update()
        self.cameraChanged.emit()

    def _showContextMenu(self, pos):
        pass

    def setGroupColor(self, groupIndex: str, color: tuple):
        self.groupColors[groupIndex] = color

    def _updateProjection(self, width: int, height: int):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        farPlane = max(self.viewSettings.maximumExtent * 2 * np.sqrt(2), 10)
        nearPlane = max(self.viewSettings.minimumExtent * 0.01, 0.01)
        gluPerspective(45.0, width / max(height, 1), nearPlane, farPlane)
        glMatrixMode(GL_MODELVIEW)

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        self._updateProjection(w, h)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._placeTimeOverlay()