from typing import Dict
import numpy as np
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GL.shaders import compileShader

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import *

from src.gui.settings import ViewSettings
from src.gui.simulator import State


class Camera:
    def __init__(self):
        self.zoom = 5.0
        self.rotationX = 45.0
        self.rotationY = 225.0
        self.minimumZoom = 1.25
        self.maximumZoom = 1000.0

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


class ObjectGroupRenderer:
    def __init__(self):
        self.vaos = {}
        self.vbos = {}
        self.colors = {}
        self.counts = {}
        self.shader = None
        self.useVaos = False

    def initialize(self):
        with open("src/assets/shaders/objects/object.vert") as f:
            vert = f.read()
        with open("src/assets/shaders/objects/object.frag") as f:
            frag = f.read()
        vertexShader = compileShader(vert, GL_VERTEX_SHADER)
        fragmentShader = compileShader(frag, GL_FRAGMENT_SHADER)
        self.shader = glCreateProgram()
        glAttachShader(self.shader, vertexShader)
        glAttachShader(self.shader, fragmentShader)
        glBindAttribLocation(self.shader, 0, "aPos")
        glLinkProgram(self.shader)
        if glGetProgramiv(self.shader, GL_LINK_STATUS) != GL_TRUE:
            log = glGetProgramInfoLog(self.shader)
            raise RuntimeError(f"ObjectGroupRenderer shader failed to link:\n{log}")
        glDeleteShader(vertexShader)
        glDeleteShader(fragmentShader)
        self.useVaos = bool(glGenVertexArrays) and bool(glBindVertexArray)

    def createGroup(self, groupIndex: str, color: tuple = (1.0, 1.0, 1.0, 0.9)):
        if groupIndex in self.vbos:
            return
        self.vbos[groupIndex] = glGenBuffers(1)
        self.colors[groupIndex] = color
        self.counts[groupIndex] = 0
        if self.useVaos:
            self.vaos[groupIndex] = glGenVertexArrays(1)
            self._configureVao(self.vaos[groupIndex], self.vbos[groupIndex])

    @staticmethod
    def _configureVao(vao, vbo):
        glBindVertexArray(vao)
        glBindBuffer(GL_ARRAY_BUFFER, vbo)
        glEnableVertexAttribArray(0)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0)

    def updateGroupPositions(self, groupIndex: str, positions: np.ndarray):
        if groupIndex not in self.vbos:
            self.createGroup(groupIndex)
        positions = np.asarray(positions, dtype=np.float32)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("Positions must be Nx3 array")
        glBindBuffer(GL_ARRAY_BUFFER, self.vbos[groupIndex])
        glBufferData(GL_ARRAY_BUFFER, positions.nbytes, positions, GL_DYNAMIC_DRAW)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        self.counts[groupIndex] = len(positions)

    def renderObjectsGroup(self, groupIndex: str, pointSize: float = 3.0):
        if groupIndex not in self.vbos or self.counts.get(groupIndex, 0) == 0:
            return
        glUseProgram(self.shader)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_PROGRAM_POINT_SIZE)
        color = self.colors[groupIndex]
        glUniform4f(glGetUniformLocation(self.shader, "uColor"), *color)
        glUniform1f(glGetUniformLocation(self.shader, "uPointSize"), pointSize)
        self._bindGroupBuffer(groupIndex)
        glDrawArrays(GL_POINTS, 0, self.counts[groupIndex])
        self._unbindGroupBuffer()
        glUseProgram(0)

    def renderAll(self, pointSize: float = 3.0):
        for groupIndex in list(self.vbos.keys()):
            self.renderObjectsGroup(groupIndex, pointSize)

    def _bindGroupBuffer(self, groupIndex):
        if self.useVaos:
            glBindVertexArray(self.vaos[groupIndex])
            return
        else:
            glBindBuffer(GL_ARRAY_BUFFER, self.vbos[groupIndex])
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)

    def _unbindGroupBuffer(self):
        if self.useVaos:
            glBindVertexArray(0)
        else:
            glDisableVertexAttribArray(0)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def cleanup(self):
        for vaoDict in [self.vaos]:
            for vao in vaoDict.values():
                glDeleteVertexArrays(1, [vao])
        for vboDict in [self.vbos]:
            for vbo in vboDict.values():
                glDeleteBuffers(1, [vbo])
        if self.shader:
            glDeleteProgram(self.shader)


class GridRenderer:
    def __init__(self):
        self.minimumExtent = 1.0
        self.maximumExtent = 1000.0
        self.linesPerHalfAxis = 10

    def initialize(self):
        pass

    def render(self, cameraZoom):
        extent = self._niceGridExtent(cameraZoom * 2.5)
        extent = max(self.minimumExtent, min(self.maximumExtent, extent))
        step = self._gridStep(extent)
        glPushMatrix()
        try:
            glUseProgram(0)
            glDisable(GL_TEXTURE_2D)
            glDisable(GL_LIGHTING)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glLineWidth(1.0)
            glBegin(GL_LINES)
            value = -extent
            while value <= extent + 1e-9:
                isAxis = abs(value) < 1e-9
                isMajor = self._isMajorGridLine(value, step)
                if isAxis:
                    glColor4f(0.85, 0.85, 0.85, 0.55)
                elif isMajor:
                    glColor4f(0.55, 0.65, 0.75, 0.22)
                else:
                    glColor4f(0.45, 0.50, 0.55, 0.10)
                glVertex3f(-extent, value, 0.0)
                glVertex3f(extent, value, 0.0)
                glVertex3f(value, -extent, 0.0)
                glVertex3f(value, extent, 0.0)
                value += step
            glEnd()
            self._drawLabels(extent, step)
        finally:
            glDisable(GL_BLEND)
            glDepthMask(GL_TRUE)
            glEnable(GL_DEPTH_TEST)
            glColor4f(1.0, 1.0, 1.0, 1.0)
            glPopMatrix()


    @staticmethod
    def _niceGridExtent(value):
        if value <= 1.0:
            return 1.0
        exponent = np.floor(np.log10(value))
        base = 10 ** exponent
        scaled = value / base
        if scaled <= 2.5:
            return base
        if scaled <= 7.5:
            return 5.0 * base
        return 10 * base

    def _gridStep(self, extent):
        return extent / self.linesPerHalfAxis

    @staticmethod
    def _isMajorGridLine(value, step):
        if abs(value) < 1e-9:
            return True
        majorStep = step * 5.0
        ratio = value / majorStep
        return abs(ratio - round(ratio)) < 1e-6

    def _drawLabels(self, extent, step):
        viewModel = (GLdouble * 16)()
        viewProjection = (GLdouble * 16)()
        viewport = (GLint * 4)()
        glGetDoublev(GL_MODELVIEW_MATRIX, viewModel)
        glGetDoublev(GL_PROJECTION_MATRIX, viewProjection)
        glGetIntegerv(GL_VIEWPORT, viewport)
        labelStep = step * 5.0
        value = -extent
        while value <= extent + 1e-9:
            if abs(value) > 1e-9 and self._isMajorGridLine(value, step):
                xLabelPosition, yLabelPosition = np.array([value, 0.0, 0.0], dtype=float), np.array([0.0, value, 0.0], dtype=float)
                xAlpha, yAlpha = 1.0, 1.0
                self._drawWorldLabel(xLabelPosition[0], xLabelPosition[1], xLabelPosition[2], self._formatGridLabel(value), viewModel, viewProjection, viewport, alpha=xAlpha, baseAlpha=0.75)
                self._drawWorldLabel(yLabelPosition[0], yLabelPosition[1], yLabelPosition[2], self._formatGridLabel(value), viewModel, viewProjection, viewport, alpha=yAlpha, baseAlpha=0.75)
            value += labelStep

    @staticmethod
    def _formatGridLabel(value):
        if abs(value) >= 10.0:
            return f"{value:.0f}"
        return f"{value:.1f}"

    @staticmethod
    def _drawWorldLabel(x, y, z, text, viewModel, viewProjection, viewport, alpha=1.0, color=(0.75, 0.80, 0.85), baseAlpha=1.0):
        xWindow, yWindow, zWindow = gluProject(x, y, z, viewModel, viewProjection, viewport)
        if zWindow <= 0.0 or zWindow >= 1.0:
            return
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, viewport[2], 0, viewport[3], -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        try:
            glUseProgram(0)
            glDisable(GL_TEXTURE_2D)
            glDisable(GL_LIGHTING)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glColor4f(color[0], color[1], color[2], baseAlpha * alpha)
            glRasterPos2f(xWindow + 4, yWindow + 4)
            for char in text:
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, ord(char))
        finally:
            glMatrixMode(GL_MODELVIEW)
            glPopMatrix()
            glMatrixMode(GL_PROJECTION)
            glPopMatrix()
            glMatrixMode(GL_MODELVIEW)


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
        for groupIndex in list(self.objectsRenderer.vbos.keys()):
            if groupIndex == 'barycenter':
                continue
            self.objectsRenderer.renderObjectsGroup(groupIndex, pointSize=4.0)
        if self.showBarycenter:
            print(self.showBarycenter)
            if "barycenter" not in self.objectsRenderer.vbos:
                self.objectsRenderer.createGroup("barycenter", (1.0, 0.15, 0.15, 1.0))
            self.objectsRenderer.updateGroupPositions("barycenter", self.massCenter.reshape(1, 3))
            self.objectsRenderer.renderObjectsGroup("barycenter", pointSize=12.0)


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
        gluPerspective(45, w / max(h, 1), 0.1, 2000 * np.sqrt(2))
        glMatrixMode(GL_MODELVIEW)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._placeTimeOverlay()