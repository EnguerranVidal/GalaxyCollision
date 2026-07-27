from dataclasses import asdict, dataclass, field
from typing import Any

from src.core.engine.parameters import SimulatorParameters


@dataclass
class WindowGeometry:
    x: int = 300
    y: int = 300
    width: int = 1200
    height: int = 600

    @classmethod
    def fromDict(cls, data=None):
        data = data or {}
        return cls(data.get("X", 300), data.get("Y", 300), data.get("WIDTH", 1200), data.get("HEIGHT", 600))

    def toDict(self):
        return {"X": self.x, "Y": self.y, "WIDTH": self.width, "HEIGHT": self.height}


@dataclass
class WindowSettings:
    maximized: bool = False
    geometry: WindowGeometry = field(default_factory=WindowGeometry)

    @classmethod
    def fromDict(cls, data=None):
        data = data or {}
        return cls(bool(data.get("MAXIMIZED", False)), WindowGeometry.fromDict(data.get("GEOMETRY")))

    def toDict(self):
        return {"MAXIMIZED": self.maximized, "GEOMETRY": self.geometry.toDict()}


@dataclass
class UiSettings:
    window: WindowSettings = field(default_factory=WindowSettings)
    parameters: SimulatorParameters = field(default_factory=SimulatorParameters)

    @classmethod
    def fromDict(cls, data: dict | None):
        data = data or {}
        return cls(window=WindowSettings.fromDict(data.get("WINDOW")),
                   parameters=SimulatorParameters.fromDict(data.get("PARAMETERS")))

    def toDict(self):
        return {"WINDOW": self.window.toDict(), "PARAMETERS": self.parameters.toDict()}
