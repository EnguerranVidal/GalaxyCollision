#include "EulerExplicitIntegrator.h"

#include "particles/ParticleGroup.h"
#include "calculators/Calculator.h"


EulerExplicitIntegrator::EulerExplicitIntegrator(float timeStep): Integrator(timeStep) {}

void EulerExplicitIntegrator::step(ParticleGroup& particles, Calculator& calculator)
{
    std::vector<Vector3> positions = particles.positions;
    std::vector<Vector3> velocities = particles.velocities;
    std::vector<Vector3> accelerations = calculator.computeAccelerations(particles);
    for (int i = 0; i < particles.nbParticles; i++)
    {
        particles.velocities[i].x = velocities[i].x + accelerations[i].x * timeStep;
        particles.velocities[i].y = velocities[i].y + accelerations[i].y * timeStep;
        particles.velocities[i].z = velocities[i].z + accelerations[i].z * timeStep;
        particles.positions[i].x = positions[i].x + velocities[i].x * timeStep;
        particles.positions[i].y = positions[i].y + velocities[i].y * timeStep;
        particles.positions[i].z = positions[i].z + velocities[i].z * timeStep;
    }


    particles.accelerations =
        accelerations;
}