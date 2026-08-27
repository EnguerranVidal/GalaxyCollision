import numpy as np
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GL import *

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import *

from src.gui.visualizers.renderers import ObjectGroupRenderer, GridRenderer
from src.gui.settings import ViewSettings
from src.gui.solver.simulator import State


class Camera:
    def __init__(self):
        self.zoom = 5.0
        self.rotationX = 45.0
        self.rotationY = 225.0
        self.minimumZoom = 1.25
        self.maximumZoom = 5000.0

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


class Universe3dViewWidget(QOpenGLWidget):
    cameraChanged = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        glutInit()
        self.setMouseTracking(True)
        self.camera = Camera()
        self.showBarycenter = False
        self.centerOnBarycenter = False
        self.massCenter = np.zeros(3, dtype=np.float32)
        self.objectSpotData, self.pendingObjectBufferUpdates = {}, {}
        self.objectsRenderer = ObjectGroupRenderer()
        self.gridRenderer = GridRenderer()
        self.pendingObjectBufferUpdates = {}
        self.groupColors = {}
        self.lastPosX, self.lastPosY = 0, 0
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._showContextMenu)
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
        self.objectsRenderer.initialize()

    def paintGL(self):
        self._uploadPendingObjectBuffers()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.camera.apply()
        self.gridRenderer.render(self.camera.zoom)
        if self.centerOnBarycenter:
            glTranslatef(-float(self.massCenter[0]), -float(self.massCenter[1]), -float(self.massCenter[2]))
        glDisable(GL_LIGHTING)
        self.objectsRenderer.renderAll(pointSize=4.0, skip={"barycenter"})
        if self.showBarycenter:
            if "barycenter" not in self.objectsRenderer.vbos:
                self.objectsRenderer.createGroup("barycenter", (1.0, 0.15, 0.15, 1.0))
            self.objectsRenderer.renderAll(pointSize=4.0, skip={"barycenter"}, refDistance=max(self.camera.zoom, 0.5), minSize=1.5, maxSize=16.0)
            glUseProgram(self.objectsRenderer.shader)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glEnable(GL_PROGRAM_POINT_SIZE)
            glEnable(GL_POINT_SPRITE)
            glDisable(GL_DEPTH_TEST)
            self.objectsRenderer.renderObjectsGroup("barycenter", pointSize=18.0)
            self.objectsRenderer.unbindGroupBuffer()
            glEnable(GL_DEPTH_TEST)
            glDisable(GL_POINT_SPRITE)
            glUseProgram(0)


    def updateState(self, state: State):
        self.pendingObjectBufferUpdates.clear()
        for groupIndex, positionArray in state.positions.items():
            if groupIndex not in self.groupColors:
                self.groupColors[groupIndex] = (0.9, 0.7, 0.3, 0.95)
            self.pendingObjectBufferUpdates[groupIndex] = positionArray
        self.massCenter = np.asarray(state.massCenter, dtype=np.float32).reshape(3)
        self.timeOverlay.setText(f"Time: {state.time:.3f}")
        self.timeOverlay.adjustSize()
        self._placeTimeOverlay()
        self.update()

    def setShowBarycenter(self, enabled: bool):
        self.showBarycenter = bool(enabled)
        if not self.showBarycenter and "barycenter" in self.objectsRenderer.counts:
            self.objectsRenderer.counts["barycenter"] = 0
        self.update()

    def setCenterOnBarycenter(self, enabled: bool):
        self.centerOnBarycenter = bool(enabled)
        self.update()

    def setViewSettings(self, viewSettings: ViewSettings):
        self.showBarycenter = bool(viewSettings.showBarycenter)
        self.centerOnBarycenter = bool(viewSettings.centerOnBarycenter)
        self.update()

    def _placeTimeOverlay(self):
        margin = 12
        x, y = self.width() - self.timeOverlay.width() - margin, self.height() - self.timeOverlay.height() - margin
        self.timeOverlay.move(max(0, x), max(0, y))
        self.timeOverlay.raise_()

    def _uploadPendingObjectBuffers(self):
        if not self.pendingObjectBufferUpdates:
            return
        for groupIndex, positions in self.pendingObjectBufferUpdates.items():
            if groupIndex not in self.objectsRenderer.colors:
                self.objectsRenderer.createGroup(groupIndex, self.groupColors.get(groupIndex, (1.0, 1.0, 1.0, 0.9)))
            self.objectsRenderer.updateGroupPositions(groupIndex, positions)
        self.pendingObjectBufferUpdates.clear()

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.lastPosX = event.x()
            self.lastPosY = event.y()

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.LeftButton:
            dx = event.x() - self.lastPosX
            dy = event.y() - self.lastPosY
            self.camera.rotate(dx, dy)
            self.lastPosX = event.x()
            self.lastPosY = event.y()
            self.update()
            self.cameraChanged.emit()

    def wheelEvent(self, event):
        delta = event.angleDelta().y() / 120.0
        self.camera.zoomBy(delta)
        self.update()
        self.cameraChanged.emit()

    def resetView(self):
        self.camera = Camera()
        self.update()
        self.cameraChanged.emit()

    def _showContextMenu(self, pos):
        pass

    def setGroupColor(self, groupIndex: str, color: tuple):
        self.groupColors[groupIndex] = color

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, w / max(h, 1), 0.1, 10000 * np.sqrt(2))
        glMatrixMode(GL_MODELVIEW)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._placeTimeOverlay()