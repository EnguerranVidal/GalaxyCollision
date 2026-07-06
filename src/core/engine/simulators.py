import cupy as cp
from src.core.engine.integrators import EulerExplicit, RK4
from src.core.engine.calculators import BarnesHutCalculator
from src.core.engine.initializers import Initializer
from src.core.engine.parameters import SimulatorParameters


class NBodySimulator:
    def __init__(self, params: SimulatorParameters):
        self.params = params
        self.timeStep = params.timeStep
        self.gravitationalConstant = params.gravitationalConstant
        self.is3D = self._getIs3D()
        self.initializer = Initializer(gravitationalConstant=self.gravitationalConstant)
        self.calculator = BarnesHutCalculator(theta=params.theta, gravitationalConstant=self.gravitationalConstant, is3D=self.is3D)
        self.integrator = self._createIntegrator(params.integratorType, self.timeStep)
        self.positions, self.velocities = None, None
        self.masses = None
        self.time = 0.0

    def _getIs3D(self):
        if self.params.distributionType.lower() == "basic":
            return self.params.basicDistributionParameters.is3D
        else:
            return self.params.galaxyDistributionParameters.is3D

    @staticmethod
    def _createIntegrator(integratorType: str, timeStep: float):
        integratorType = integratorType.upper()
        if integratorType == "EULER_EXPLICIT":
            return EulerExplicit(timeStep)
        else:
            return RK4(timeStep)

    def initialize(self):
        distributionType = self.params.distributionType.lower()
        if distributionType == "basic":
            p = self.params.basicDistributionParameters
            self.positions, self.velocities, self.masses = self.initializer.basicDistribution(numParticles=p.numParticles, positionScale=p.positionScale, velocityScale=p.velocityScale, is3D=p.is3D)
        elif distributionType == "galaxy":
            p = self.params.galaxyDistributionParameters
            self.positions, self.velocities, self.masses = self.initializer.galaxyDistribution(numParticles=p.numParticles, totalMass=p.totalMass, radius=p.radius, height=p.height, is3D=p.is3D)
        else:
            raise ValueError("Unsupported distributionType")
        print(f"Initialized {distributionType} simulation with {self.params.numParticles if hasattr(self.params, 'numParticles') else 'N'} particles.")

    def step(self):
        if self.positions is None:
            raise ValueError("Particles not initialized. Call initializeBasic() or initializeGalaxy() first.")
        self.positions, self.velocities = self.integrator.step(self.positions, self.velocities, self.masses, self.calculator)
        self.time += self.timeStep
        return self.positions, self.velocities

    def runSimulation(self, numSteps: int = 1000, printEvery: int = 100):
        print("Simulation started!")
        for step in range(numSteps):
            self.step()
            if (step + 1) % printEvery == 0:
                print(f"Step {step + 1}/{numSteps} | Time: {self.time:.4f}")
        print("Simulation completed!")
        return self.positions, self.velocities, self.time
