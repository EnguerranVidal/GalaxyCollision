import numpy as np
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader

import cupy as cp
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import *


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
        rotationX = cp.deg2rad(self.rotationX)
        rotationY = cp.deg2rad(self.rotationY)
        rotationFix = cp.deg2rad(-90)
        camera = cp.array([0, 0, self.zoom], dtype=float)
        xRotation = cp.array([[1, 0, 0], [0, cp.cos(-rotationX), -cp.sin(-rotationX)], [0, cp.sin(-rotationX), cp.cos(-rotationX)]])
        yRotation = cp.array([[cp.cos(-rotationY), 0, cp.sin(-rotationY)], [0, 1, 0], [-cp.sin(-rotationY), 0, cp.cos(-rotationY)]])
        fixRotation = cp.array([[1, 0, 0], [0, cp.cos(-rotationFix), -cp.sin(-rotationFix)], [0, cp.sin(-rotationFix), cp.cos(-rotationFix)]])
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

    def updateGroupPositions(self, groupIndex: str, positions: cp.ndarray):
        if groupIndex not in self.vbos:
            self.createGroup(groupIndex)
        positions = cp.asarray(positions, dtype=cp.float32)
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


class Universe3dViewWidget(QOpenGLWidget):
    cameraChanged = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        glutInit()
        self.setMouseTracking(True)
        self.camera = Camera()
        self.objectSpotData, self.pendingObjectBufferUpdates = {}, {}
        self.objectsRenderer = ObjectGroupRenderer()
        self.pendingObjectBufferUpdates = {}
        self.groupColors = {}
        self.lastPosX, self.lastPosY = 0, 0
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._showContextMenu)

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
        glDisable(GL_LIGHTING)
        self.objectsRenderer.renderAll(pointSize=4.0)

    def updateData(self, positions: dict):
        self.pendingObjectBufferUpdates.clear()
        for groupIndex, positionArray in positions.items():
            if groupIndex not in self.groupColors:
                self.groupColors[groupIndex] = (0.9, 0.7, 0.3, 0.95)
            self.pendingObjectBufferUpdates[groupIndex] = positionArray
        self.update()

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