import os
import json
from typing import Optional

from PyQt5.QtWidgets import *

from src.gui.settings import UiSettings, WindowGeometry


class MainWindow(QMainWindow):
    def __init__(self, currentDir: str):
        super().__init__()
        self.icons = {}
        self.settings: Optional[UiSettings] = None
        # FOLDER PATHS & SETTINGS
        self.currentDir = currentDir
        self.settingsPath = os.path.join(self.currentDir, 'settings.json')
        self._checkEnvironment()
        self.loadSettings()
        self._restoreWindow()

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
            windowGeometry = self.settings.window.geometry
            self.setGeometry(windowGeometry.x, windowGeometry.y, windowGeometry.width, windowGeometry.height)


    def closeEvent(self, event):
        self.settings.window.maximized = self.isMaximized()
        if not self.isMaximized():
            g = self.geometry()
            self.settings.window.geometry = WindowGeometry.fromDict({'X': g.x(), 'Y': g.y(), 'WIDTH': g.width(), 'HEIGHT': g.height()})
        self.saveSettings()
        event.accept()