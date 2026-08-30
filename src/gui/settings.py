from dataclasses import asdict, dataclass, field

from src.gui.solver.parameters import SimulatorParameters


@dataclass
class WindowSettings:
    maximized: bool = False
    x: int = 300
    y: int = 300
    width: int = 1200
    height: int = 600

    @classmethod
    def fromDict(cls, data=None):
        data = data or {}
        return cls(bool(data.get("MAXIMIZED", False)), data.get("X", 300), data.get("Y", 300), data.get("WIDTH", 1200), data.get("HEIGHT", 600))

    def toDict(self):
        return {"MAXIMIZED": self.maximized, "X": self.x, "Y": self.y, "WIDTH": self.width, "HEIGHT": self.height}


@dataclass
class ViewSettings:
    showBarycenter: bool = False
    centerOnBarycenter: bool = False
    showVelocityVectors: bool = False
    velocityVectorLength: float = 0.25
    referenceVelocity: float = 1E1
    showAccelerationVectors: bool = False
    accelerationVectorLength: float = 0.25
    referenceAcceleration: float = 1E1

    @classmethod
    def fromDict(cls, data=None):
        data = data or {}
        return cls(
            showBarycenter=bool(data.get("SHOW_BARYCENTER", False)),
            centerOnBarycenter=bool(data.get("CENTER_ON_BARYCENTER", False)),
            showVelocityVectors=bool(data.get("SHOW_VELOCITY_VECTORS", False)),
            velocityVectorLength=float(data.get("VELOCITY_VECTOR_LENGTH", 0.25)),
            referenceVelocity=float(data.get("REFERENCE_VELOCITY", 1E1)),
            showAccelerationVectors=bool(data.get("SHOW_ACCELERATION_VECTORS", False)),
            accelerationVectorLength=float(data.get("ACCELERATION_VECTOR_LENGTH", 0.25)),
            referenceAcceleration=float(data.get("REFERENCE_ACCELERATION", 1E1)),
        )

    def toDict(self):
        return {
            "SHOW_BARYCENTER": self.showBarycenter,
            "CENTER_ON_BARYCENTER": self.centerOnBarycenter,
            "SHOW_VELOCITY_VECTORS": self.showVelocityVectors,
            "VELOCITY_VECTOR_LENGTH": self.velocityVectorLength,
            "REFERENCE_VELOCITY": self.referenceVelocity,
            "SHOW_ACCELERATION_VECTORS": self.showAccelerationVectors,
            "ACCELERATION_VECTOR_LENGTH": self.accelerationVectorLength,
            "REFERENCE_ACCELERATION": self.referenceAcceleration,
        }


@dataclass
class UiSettings:
    window: WindowSettings = field(default_factory=WindowSettings)
    parameters: SimulatorParameters = field(default_factory=SimulatorParameters)
    view: ViewSettings = field(default_factory=ViewSettings)

    @classmethod
    def fromDict(cls, data: dict | None):
        data = data or {}
        return cls(window=WindowSettings.fromDict(data.get("WINDOW", WindowSettings())),
                   parameters=SimulatorParameters.fromDict(data.get("PARAMETERS", SimulatorParameters())),
                   view=ViewSettings.fromDict(data.get("VIEW", ViewSettings())))

    def toDict(self):
        return {"WINDOW": self.window.toDict(), "PARAMETERS": self.parameters.toDict(), "VIEW": self.view.toDict()}
