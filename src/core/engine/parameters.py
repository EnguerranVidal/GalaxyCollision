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


@dataclass
class GalaxyDistributionParameters:
    numParticles: int = 5000
    totalMass: float = 1000.0
    radius: float = 10.0
    height: float = 2.0
    is3D: bool = True


@dataclass
class SimulatorParameters:
    timeStep: float = 0.001
    theta: float = 0.5
    is3D = True
    gravitationalConstant: float = 1.0
    integratorType: str = "RK4"
    distributionType: str = "GALAXY"
    calculatorType: str = "BARNES-HUT"
    basicDistributionParameters: Optional[BasicDistributionParameters] = None
    galaxyDistributionParameters: Optional[GalaxyDistributionParameters] = None

    def __post_init__(self):
        if self.distributionType.lower() == "basic" and self.basicDistributionParameters is None:
            self.basicDistributionParameters = BasicDistributionParameters(is3D=self.is3D)
        elif self.distributionType.lower() == "galaxy" and self.galaxyDistributionParameters is None:
            self.galaxyDistributionParameters = GalaxyDistributionParameters(is3D=self.is3D)
