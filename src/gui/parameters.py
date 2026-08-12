from dataclasses import dataclass, field
from typing import Optional
import engine


@dataclass
class BasicDistributionParameters:
    nbParticles: int = 1000
    positionScale: float = 10.0
    velocityScale: float = 1.0
    massMinimum: float = 0.1
    massMaximum: float = 1.0

    @classmethod
    def fromDict(cls, data: dict):
        return cls(**data)

    def toDict(self):
        return {
            "nbParticles": self.nbParticles,
            "positionScale": self.positionScale,
            "velocityScale": self.velocityScale,
            "massMinimum": self.massMinimum,
            "massMaximum": self.massMaximum,
        }

    def __eq__(self, other) -> bool:
        if not isinstance(other, BasicDistributionParameters):
            return NotImplemented
        return (
                self.nbParticles == other.nbParticles and
                abs(self.positionScale - other.positionScale) < 1e-9 and
                abs(self.velocityScale - other.velocityScale) < 1e-9 and
                abs(self.massMinimum - other.massMinimum) < 1e-9 and
                abs(self.massMaximum - other.massMaximum) < 1e-9
        )

    def __repr__(self):
        return f"BasicDistributionParameters(numParticles={self.nbParticles}, positionScale={self.positionScale})"

    def toCpp(self):
        cpp = engine.BasicDistributionParameters()
        cpp.nbParticles = self.nbParticles
        cpp.positionScale = self.positionScale
        cpp.velocityScale = self.velocityScale
        cpp.massMinimum = self.massMinimum
        cpp.massMaximum = self.massMaximum
        return cpp

    @classmethod
    def fromCpp(cls, cpp):
        return cls(
            nbParticles=cpp.nbParticles,
            positionScale=cpp.positionScale,
            velocityScale=cpp.velocityScale,
            massMinimum=cpp.massMinimum,
            massMaximum=cpp.massMaximum,
        )


@dataclass
class GalaxyDistributionParameters:
    nbParticles: int = 5000
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
            "nbParticles": self.nbParticles,
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
                self.nbParticles == other.nbParticles and
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
        return f"GalaxyDistributionParameters(numParticles={self.nbParticles}, radius={self.radius}, height={self.height})"

    def toCpp(self):
        cpp = engine.GalaxyDistributionParameters()
        cpp.nbParticles = self.nbParticles
        cpp.totalMass = self.totalMass
        cpp.radius = self.radius
        cpp.height = self.height
        cpp.bulgeFraction = self.bulgeFraction
        cpp.haloFraction = self.haloFraction
        cpp.velocityDispersion = self.velocityDispersion
        cpp.plummerRadius = self.plummerRadius
        cpp.haloRadius = self.haloRadius
        return cpp

    @classmethod
    def fromCpp(cls, cpp):
        return cls(
            nbParticles=cpp.nbParticles,
            totalMass=cpp.totalMass,
            radius=cpp.radius,
            height=cpp.height,
            bulgeFraction=cpp.bulgeFraction,
            haloFraction=cpp.haloFraction,
            velocityDispersion=cpp.velocityDispersion,
            plummerRadius=cpp.plummerRadius,
            haloRadius=cpp.haloRadius,
        )


@dataclass
class SimulatorParameters:
    name: str = "Untitled Simulation"
    timeStep: float = 0.001
    theta: float = 0.5
    seed: int = 0
    device: str = "CPU"
    gravitationalConstant: float = 1.0
    integratorType: str = "RK4"
    distributionType: str = "GALAXY"
    calculatorType: str = "BARNES_HUT"
    endless: bool = True
    maxTime: float = 1000.0
    saveResults: bool = False
    basicDistributionParameters: Optional[BasicDistributionParameters] = field(default_factory=BasicDistributionParameters)
    galaxyDistributionParameters: Optional[GalaxyDistributionParameters] = field(default_factory=GalaxyDistributionParameters)

    @classmethod
    def fromDict(cls, data: dict):
        if data is None:
            return cls()
        data = dict(data)
        basic = data.get("basicDistributionParameters")
        if basic is not None:
            data["basicDistributionParameters"] = basic if isinstance(basic, BasicDistributionParameters) else BasicDistributionParameters.fromDict(basic)
        galaxy = data.get("galaxyDistributionParameters")
        if galaxy is not None:
            data["galaxyDistributionParameters"] = galaxy if isinstance(galaxy, GalaxyDistributionParameters) else GalaxyDistributionParameters.fromDict(galaxy)
        return cls(**data)

    def toDict(self):
        return {
            "name": self.name,
            "timeStep": self.timeStep,
            "theta": self.theta,
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
        if not isinstance(other, SimulatorParameters):
            return NotImplemented
        mainEquality = (
                self.name == other.name and
                abs(self.timeStep - other.timeStep) < 1e-9 and
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
        distributionType = self.galaxyDistributionParameters if self.distributionType.upper() == "GALAXY" else self.basicDistributionParameters
        nbParticles = distributionType.nbParticles if distributionType else 0
        return f"SimulatorParameters(name='{self.name}', seed={self.seed}, device={self.device}, distribution={self.distributionType}, particles={nbParticles})"


    def toCpp(self):
        cpp = engine.SimulatorParameters()
        cpp.name = self.name
        cpp.timeStep = self.timeStep
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

    @classmethod
    def fromCpp(cls, cpp):
        return cls(
            name = cpp.name,
            timeStep = cpp.timeStep,
            theta = cpp.theta,
            seed = cpp.seed,
            device = cpp.device,
            gravitationalConstant = cpp.gravitationalConstant,
            integratorType = cpp.integratorType,
            distributionType = cpp.distributionType,
            calculatorType = cpp.calculatorType,
            endless = cpp.endless,
            maxTime = cpp.maxTime,
            saveResults = cpp.saveResults,
            galaxyDistributionParameters = cpp.galaxyDistributionParameters.fromCpp(),
            basicDistributionParameters = cpp.basicDistributionParameters.fromCpp()
        )