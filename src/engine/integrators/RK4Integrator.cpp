#include "RK4Integrator.h"

#include "math/Vector3.h"
#include "particles/ParticleGroup.h"
#include "calculator/Calculator.h"


namespace
{

Vector3 weightedAverage(const Vector3& k1, const Vector3& k2, const Vector3& k3, const Vector3& k4)
{
    return {
        (k1.x + 2.0f * k2.x + 2.0f * k3.x + k4.x) / 6.0f,
        (k1.y + 2.0f * k2.y + 2.0f * k3.y + k4.y) / 6.0f,
        (k1.z + 2.0f * k2.z + 2.0f * k3.z + k4.z) / 6.0f
    };
}
}

RK4::RK4(float timeStep): Integrator(timeStep) {}

void RK4::step(ParticleGroup& particles, Calculator& calculator)
{
    std::vector<Vector3> positions = particles.positions;
    std::vector<Vector3> velocities = particles.velocities;

    // K1 COMPUTATIONS
    std::vector<Vector3> k1Accelerations = calculator.computeAccelerations(particles);
    std::vector<Vector3> k1Velocities(particles.nbParticles);
    std::vector<Vector3> k1Positions(particles.nbParticles);
    for (int i = 0; i < particles.nbParticles; i++)
    {
        k1Velocities[i] = k1Accelerations[i] * timeStep;
        k1Positions[i] = velocities[i] * timeStep;
    }

    // K2 COMPUTATIONS
    for (int i = 0; i < particles.nbParticles; i++)
    {
        particles.positions[i] = positions[i] + k1Positions[i] * 0.5f;
        particles.velocities[i] = velocities[i] + k1Velocities[i] * 0.5f;
    }
    std::vector<Vector3> k2Accelerations = calculator.computeAccelerations(particles);
    std::vector<Vector3> k2Velocities(particles.nbParticles);
    std::vector<Vector3> k2Positions(particles.nbParticles);
    for (int i = 0; i < particles.nbParticles; i++)
    {
        k2Velocities[i] = k2Accelerations[i] * timeStep;
        k2Positions[i] = (velocities[i] + k1Velocities[i] * 0.5f) * timeStep;
    }

    // K3 COMPUTATIONS
    for (int i = 0; i < particles.nbParticles; i++)
    {
        particles.positions[i] = positions[i] + k2Positions[i] * 0.5f;
        particles.velocities[i] = velocities[i] + k2Velocities[i] * 0.5f;
    }
    std::vector<Vector3> k3Accelerations = calculator.computeAccelerations(particles);
    std::vector<Vector3> k3Velocities(particles.nbParticles);
    std::vector<Vector3> k3Positions(particles.nbParticles);
    for (int i = 0; i < particles.nbParticles; i++)
    {
        k3Velocities[i] = k3Accelerations[i] * timeStep;
        k3Positions[i] = (velocities[i] + k2Velocities[i] * 0.5f) * timeStep;
    }

    // K4 COMPUTATIONS
    for (int i = 0; i < particles.nbParticles; i++)
    {
        particles.positions[i] = positions[i] + k3Positions[i];
        particles.velocities[i] = velocities[i] + k3Velocities[i];
    }
    std::vector<Vector3> k4Accelerations = calculator.computeAccelerations(particles);
    std::vector<Vector3> k4Velocities(particles.nbParticles);
    std::vector<Vector3> k4Positions(particles.nbParticles);
    for (int i = 0; i < particles.nbParticles; i++)
    {
        k4Velocities[i] = k4Accelerations[i] * timeStep;
        k4Positions[i] = (velocities[i] + k3Velocities[i]) * timeStep;
    }
    for (int i = 0; i < particles.nbParticles; i++)
    {
        particles.positions[i] = positions[i] + weightedAverage(k1Positions[i], k2Positions[i], k3Positions[i], k4Positions[i]);
        particles.velocities[i] = velocities[i] + weightedAverage(k1Velocities[i], k2Velocities[i], k3Velocities[i], k4Velocities[i]);
    }
    particles.accelerations = k4Accelerations;
}