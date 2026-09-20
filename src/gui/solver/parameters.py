from dataclasses import dataclass, field
from typing import Optional
import engine


@dataclass
class BasicDistributionParameters:
    positionScale: float = 10.0
    velocityScale: float = 1.0
    massMinimum: float = 0.1
    massMaximum: float = 1.0

    @classmethod
    def fromDict(cls, data: dict):
        return cls(**data)

    def toDict(self):
        return {
            "positionScale": self.positionScale,
            "velocityScale": self.velocityScale,
            "massMinimum": self.massMinimum,
            "massMaximum": self.massMaximum,
        }

    def __eq__(self, other) -> bool:
        if not isinstance(other, BasicDistributionParameters):
            return NotImplemented
        return (
                abs(self.positionScale - other.positionScale) < 1e-9 and
                abs(self.velocityScale - other.velocityScale) < 1e-9 and
                abs(self.massMinimum - other.massMinimum) < 1e-9 and
                abs(self.massMaximum - other.massMaximum) < 1e-9
        )

    def __repr__(self):
        return f"BasicDistributionParameters(positionScale={self.positionScale})"

    def toCpp(self):
        cpp = engine.BasicDistributionParameters()
        cpp.positionScale = self.positionScale
        cpp.velocityScale = self.velocityScale
        cpp.massMinimum = self.massMinimum
        cpp.massMaximum = self.massMaximum
        return cpp


@dataclass
class GalaxyDistributionParameters:
    totalMass: float = 1000.0
    radius: float = 15.0
    height: float = 2.0
    bulgeFraction: float = 0.2
    haloFraction: float = 0.2
    velocityDispersion: float = 0.1
    plummerRadius: float = 3.0
    haloRadius: float = 5.0

    @classmethod
    def fromDict(cls, data: dict):
        return cls(**data)

    def toDict(self):
        return {
            "totalMass": self.totalMass,
            "radius": self.radius,
            "height": self.height,
            "bulgeFraction": self.bulgeFraction,
            "haloFraction": self.haloFraction,
            "velocityDispersion": self.velocityDispersion,
            "plummerRadius": self.plummerRadius,
            "haloRadius": self.haloRadius,
        }

    def __eq__(self, other) -> bool:
        if not isinstance(other, GalaxyDistributionParameters):
            return NotImplemented
        return (
                abs(self.totalMass - other.totalMass) < 1e-9 and
                abs(self.radius - other.radius) < 1e-9 and
                abs(self.height - other.height) < 1e-9 and
                abs(self.bulgeFraction - other.bulgeFraction) < 1e-9 and
                abs(self.haloFraction - other.haloFraction) < 1e-9 and
                abs(self.velocityDispersion - other.velocityDispersion) < 1e-9 and
                abs(self.plummerRadius - other.plummerRadius) < 1e-9 and
                abs(self.haloRadius - other.haloRadius) < 1e-9
        )

    def __repr__(self):
        return f"GalaxyDistributionParameters(radius={self.radius}, height={self.height})"

    def toCpp(self):
        cpp = engine.GalaxyDistributionParameters()
        cpp.totalMass = self.totalMass
        cpp.radius = self.radius
        cpp.height = self.height
        cpp.bulgeFraction = self.bulgeFraction
        cpp.haloFraction = self.haloFraction
        cpp.velocityDispersion = self.velocityDispersion
        cpp.plummerRadius = self.plummerRadius
        cpp.haloRadius = self.haloRadius
        return cpp


@dataclass
class SolverParameters:
    name: str = "Untitled Simulation"
    timeStep: float = 5E-3
    maxFps: float = 60.0
    simulationSpeed: float = 1.0
    nbParticles: int = 10000
    theta: float = 0.5
    tileSize: int = 32
    blockSize: int = 256
    seed: int = 0
    device: str = "GPU"
    gravitationalConstant: float = 1.0
    integratorType: str = "EULER_EXPLICIT"
    distributionType: str = "GALAXY"
    calculatorType: str = "NEWTON"
    endless: bool = True
    maxTime: float = 1000.0
    saveResults: bool = False
    basicDistributionParameters: Optional[BasicDistributionParameters] = field(default_factory=BasicDistributionParameters)
    galaxyDistributionParameters: Optional[GalaxyDistributionParameters] = field(default_factory=GalaxyDistributionParameters)

    @classmethod
    def fromDict(cls, data: dict):
        data = dict(data)
        basic = data.pop("basicDistributionParameters", None)
        if isinstance(basic, dict):
            data["basicDistributionParameters"] = basic if isinstance(basic, BasicDistributionParameters) else BasicDistributionParameters.fromDict(basic)
        galaxy = data.pop("galaxyDistributionParameters", None)
        if isinstance(galaxy, dict):
            data["galaxyDistributionParameters"] = galaxy if isinstance(galaxy, GalaxyDistributionParameters) else GalaxyDistributionParameters.fromDict(galaxy)
        return cls(**data)

    def toDict(self):
        return {
            "name": self.name,
            "timeStep": self.timeStep,
            "maxFps": self.maxFps,
            "simulationSpeed": self.simulationSpeed,
            "nbParticles": self.nbParticles,
            "theta": self.theta,
            "tileSize": self.tileSize,
            "blockSize": self.blockSize,
            "seed": self.seed,
            "device": self.device,
            "gravitationalConstant": self.gravitationalConstant,
            "integratorType": self.integratorType,
            "distributionType": self.distributionType,
            "calculatorType": self.calculatorType,
            "endless": self.endless,
            "maxTime": self.maxTime,
            "saveResults": self.saveResults,
            "basicDistributionParameters": self.basicDistributionParameters.toDict() if self.basicDistributionParameters is not None else None,
            "galaxyDistributionParameters": self.galaxyDistributionParameters.toDict() if self.galaxyDistributionParameters is not None else None,
        }

    def __eq__(self, other) -> bool:
        if not isinstance(other, SolverParameters):
            return NotImplemented
        mainEquality = (
                self.name == other.name and
                abs(self.timeStep - other.timeStep) < 1e-9 and
                abs(self.maxFps - other.maxFps) < 1e-9 and
                abs(self.simulationSpeed - other.simulationSpeed) < 1e-9 and
                self.nbParticles == other.nbParticles and
                abs(self.theta - other.theta) < 1e-9 and
                self.seed == other.seed and
                self.device == other.device and
                abs(self.gravitationalConstant - other.gravitationalConstant) < 1e-9 and
                self.integratorType.upper() == other.integratorType.upper() and
                self.distributionType.upper() == other.distributionType.upper() and
                self.calculatorType.upper() == other.calculatorType.upper() and
                self.endless == other.endless and
                abs(self.maxTime - other.maxTime) < 1e-6 and
                self.saveResults == other.saveResults
        )
        if not mainEquality:
            return False
        if self.distributionType.upper() == "GALAXY":
            return self.galaxyDistributionParameters == other.galaxyDistributionParameters
        else:
            return self.basicDistributionParameters == other.basicDistributionParameters

    def __repr__(self):
        return f"SolverParameters(name='{self.name}', seed={self.seed}, device={self.device}, distribution={self.distributionType}, particles={self.nbParticles})"

    def toCpp(self):
        cpp = engine.SolverParameters()
        cpp.name = self.name
        cpp.timeStep = self.timeStep
        cpp.nbParticles = self.nbParticles
        cpp.theta = self.theta
        cpp.seed = self.seed
        cpp.device = self.device.lower()
        cpp.gravitationalConstant = self.gravitationalConstant
        cpp.integratorType = self.integratorType
        cpp.distributionType = self.distributionType
        cpp.calculatorType = self.calculatorType
        cpp.endless = self.endless
        cpp.maxTime = self.maxTime
        cpp.saveResults = self.saveResults
        if self.distributionType.upper() == "GALAXY":
            cpp.galaxyDistributionParameters = self.galaxyDistributionParameters.toCpp()
        else:
            cpp.basicDistributionParameters = self.basicDistributionParameters.toCpp()
        return cpp
