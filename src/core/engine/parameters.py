from dataclasses import dataclass
from typing import Optional


@dataclass
class BasicDistributionParameters:
    numParticles: int = 1000
    positionScale: float = 10.0
    velocityScale: float = 1.0
    massMin: float = 0.1
    massMax: float = 1.0
    is3D: bool = True

    def __eq__(self, other) -> bool:
        if not isinstance(other, BasicDistributionParameters):
            return NotImplemented
        return (
                self.numParticles == other.numParticles and
                abs(self.positionScale - other.positionScale) < 1e-9 and
                abs(self.velocityScale - other.velocityScale) < 1e-9 and
                abs(self.massMin - other.massMin) < 1e-9 and
                abs(self.massMax - other.massMax) < 1e-9 and
                self.is3D == other.is3D
        )

    def __repr__(self):
        return f"BasicDistributionParameters(numParticles={self.numParticles}, positionScale={self.positionScale}, is3D={self.is3D})"


@dataclass
class GalaxyDistributionParameters:
    numParticles: int = 5000
    totalMass: float = 1000.0
    radius: float = 10.0
    height: float = 2.0
    is3D: bool = True

    def __eq__(self, other) -> bool:
        if not isinstance(other, GalaxyDistributionParameters):
            return NotImplemented
        return (
                self.numParticles == other.numParticles and
                abs(self.totalMass - other.totalMass) < 1e-9 and
                abs(self.radius - other.radius) < 1e-9 and
                abs(self.height - other.height) < 1e-9 and
                self.is3D == other.is3D
        )

    def __repr__(self):
        return f"GalaxyDistributionParameters(numParticles={self.numParticles}, radius={self.radius}, height={self.height}, is3D={self.is3D})"


@dataclass
class SimulatorParameters:
    name: str = "Untitled Simulation"
    timeStep: float = 0.001
    theta: float = 0.5
    is3D: bool = True
    gravitationalConstant: float = 1.0
    integratorType: str = "RK4"
    distributionType: str = "GALAXY"
    calculatorType: str = "BARNES-HUT"
    endless: bool = True
    maxTime: float = 1000.0
    saveResults: bool = False
    basicDistributionParameters: Optional[BasicDistributionParameters] = None
    galaxyDistributionParameters: Optional[GalaxyDistributionParameters] = None

    def __post_init__(self):
        if self.distributionType.lower() == "basic" and self.basicDistributionParameters is None:
            self.basicDistributionParameters = BasicDistributionParameters(is3D=self.is3D)
        elif self.distributionType.lower() == "galaxy" and self.galaxyDistributionParameters is None:
            self.galaxyDistributionParameters = GalaxyDistributionParameters(is3D=self.is3D)

    def __eq__(self, other):
        if not isinstance(other, SimulatorParameters):
            return NotImplemented
        mainEquality = (
                self.name == other.name and
                abs(self.timeStep - other.timeStep) < 1e-9 and
                abs(self.theta - other.theta) < 1e-9 and
                abs(self.gravitationalConstant - other.gravitationalConstant) < 1e-9 and
                self.is3D == other.is3D and
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
        nbParticles = distributionType.numParticles if distributionType else 0
        return f"SimulatorParameters(name='{self.name}', distribution={self.distributionType}, particles={nbParticles})"
