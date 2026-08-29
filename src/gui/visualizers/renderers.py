import numpy as np
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GL.shaders import compileShader


class ParticlesRenderer:
    def __init__(self):
        self.vaos = {}
        self.vbos = {}
        self.colors = {}
        self.counts = {}
        self.capacities = {}
        self.shader = None
        self.useVaos = False
        self.locColor = -1
        self.locPointSize = -1
        self.locRefDistance = -1
        self.locMinSize = -1
        self.locMaxSize = -1

    def initialize(self):
        with open("src/assets/shaders/particles/particles.vert") as f:
            vert = f.read()
        with open("src/assets/shaders/particles/particles.frag") as f:
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
            raise RuntimeError(f"ParticlesRenderer shader failed to link:\n{log}")
        glDeleteShader(vertexShader)
        glDeleteShader(fragmentShader)
        self.locColor = glGetUniformLocation(self.shader, "uColor")
        self.locPointSize = glGetUniformLocation(self.shader, "uPointSize")
        self.locRefDistance = glGetUniformLocation(self.shader, "uRefDistance")
        self.locMinSize = glGetUniformLocation(self.shader, "uMinSize")
        self.locMaxSize = glGetUniformLocation(self.shader, "uMaxSize")
        self.useVaos = bool(glGenVertexArrays) and bool(glBindVertexArray)

    def createGroup(self, groupIndex: str, color: tuple = (1.0, 1.0, 1.0, 0.9)):
        if groupIndex in self.vbos:
            return
        self.vbos[groupIndex] = glGenBuffers(1)
        self.colors[groupIndex] = color
        self.counts[groupIndex] = 0
        self.capacities[groupIndex] = 0
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
        positions = np.ascontiguousarray(positions, dtype=np.float32)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("Positions must be Nx3 array")
        nbParticles = positions.shape[0]
        nbBytes = int(positions.nbytes)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbos[groupIndex])
        capacity = self.capacities.get(groupIndex, 0)
        if nbBytes > capacity:
            newCapacity = max(nbBytes, int(capacity * 1.5) if capacity else nbBytes)
            glBufferData(GL_ARRAY_BUFFER, newCapacity, None, GL_DYNAMIC_DRAW)
            self.capacities[groupIndex] = newCapacity
        else:
            glBufferData(GL_ARRAY_BUFFER, capacity, None, GL_DYNAMIC_DRAW)
        glBufferSubData(GL_ARRAY_BUFFER, 0, nbBytes, positions)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        self.counts[groupIndex] = nbParticles

    def renderObjectsGroup(self, groupIndex: str, pointSize: float = 3.0):
        if groupIndex not in self.vbos or self.counts.get(groupIndex, 0) == 0:
            return
        glUniform4f(self.locColor, *self.colors[groupIndex])
        glUniform1f(self.locPointSize, pointSize)
        self._bindGroupBuffer(groupIndex)
        glDrawArrays(GL_POINTS, 0, self.counts[groupIndex])

    def renderAll(self, pointSize: float = 3.0, skip=None, refDistance: float = None, minSize: float = 1.0, maxSize: float = 12.0):
        skip = skip or set()
        if refDistance is None:
            refDistance = pointSize
        glUseProgram(self.shader)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_PROGRAM_POINT_SIZE)
        glEnable(GL_POINT_SPRITE)
        glUniform1f(self.locPointSize, pointSize)
        glUniform1f(self.locRefDistance, float(refDistance))
        glUniform1f(self.locMinSize, minSize)
        glUniform1f(self.locMaxSize, maxSize)
        for groupIndex, count in self.counts.items():
            if count == 0 or groupIndex in skip:
                continue
            glUniform4f(self.locColor, *self.colors[groupIndex])
            glUniform1f(self.locPointSize, pointSize)
            self._bindGroupBuffer(groupIndex)
            glDrawArrays(GL_POINTS, 0, count)
        self.unbindGroupBuffer()
        glDisable(GL_POINT_SPRITE)
        glUseProgram(0)

    def _bindGroupBuffer(self, groupIndex):
        if self.useVaos:
            glBindVertexArray(self.vaos[groupIndex])
        else:
            glBindBuffer(GL_ARRAY_BUFFER, self.vbos[groupIndex])
            glEnableVertexAttribArray(0)
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)

    def unbindGroupBuffer(self):
        if self.useVaos:
            glBindVertexArray(0)
        else:
            glDisableVertexAttribArray(0)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def cleanup(self):
        for vao in self.vaos.values():
            glDeleteVertexArrays(1, [vao])
        for vbo in self.vbos.values():
            glDeleteBuffers(1, [vbo])
        self.vaos.clear()
        self.vbos.clear()
        self.colors.clear()
        self.counts.clear()
        self.capacities.clear()
        if self.shader:
            glDeleteProgram(self.shader)
            self.shader = None


class BarycenterRenderer:
    def __init__(self):
        self.vbo = None
        self.shader = None
        self.locColor = -1
        self.locPointSize = -1
        self.ready = False

    def initialize(self):
        self.vbo = glGenBuffers(1)
        with open("src/assets/shaders/barycenter/barycenter.vert") as f:
            vert = f.read()
        with open("src/assets/shaders/barycenter/barycenter.frag") as f:
            frag = f.read()
        vertexShader = compileShader(vert, GL_VERTEX_SHADER)
        fragmentShader = compileShader(frag, GL_FRAGMENT_SHADER)
        self.shader = glCreateProgram()
        glAttachShader(self.shader, vertexShader)
        glAttachShader(self.shader, fragmentShader)
        glBindAttribLocation(self.shader, 0, "aPos")
        glLinkProgram(self.shader)
        if glGetProgramiv(self.shader, GL_LINK_STATUS) != GL_TRUE:
            raise RuntimeError(glGetProgramInfoLog(self.shader))
        glDeleteShader(vertexShader)
        glDeleteShader(fragmentShader)
        self.locColor = glGetUniformLocation(self.shader, "uColor")
        self.locPointSize = glGetUniformLocation(self.shader, "uPointSize")
        self.ready = True

    def render(self, position: np.ndarray, pointSize: float = 14.0, color=(1.0, 0.15, 0.15, 1.0)):
        if not self.ready:
            return
        pos = np.ascontiguousarray(position, dtype=np.float32).reshape(1, 3)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, pos.nbytes, pos, GL_DYNAMIC_DRAW)
        glUseProgram(self.shader)
        glEnable(GL_PROGRAM_POINT_SIZE)
        glEnable(GL_POINT_SPRITE)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glDisable(GL_DEPTH_TEST)
        glUniform4f(self.locColor, *color)
        glUniform1f(self.locPointSize, float(pointSize))
        glEnableVertexAttribArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glDrawArrays(GL_POINTS, 0, 1)
        glDisableVertexAttribArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_POINT_SPRITE)
        glUseProgram(0)


class GridRenderer:
    def __init__(self):
        self.minimumExtent = 1.0
        self.maximumExtent = 5000.0
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
