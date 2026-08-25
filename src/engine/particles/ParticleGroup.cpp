#include "ParticleGroup.h"

#include <stdexcept>


ParticleGroup::ParticleGroup(int nbParticles, const std::string& device):
    device(device),
    nbParticles(nbParticles)
{
    if (device != "cpu" && device != "gpu") {throw std::runtime_error("Unknown device: " + device);}
    positions.resize(nbParticles);
    velocities.resize(nbParticles);
    accelerations.resize(nbParticles);
    masses.resize(nbParticles, 1.0f);
}

void ParticleGroup::setPositions(const std::vector<Vector3>& newPositions)
{
    positions = newPositions;
    nbParticles = static_cast<int>(positions.size());
}

void ParticleGroup::setVelocities(const std::vector<Vector3>& newVelocities) {velocities = newVelocities;}

void ParticleGroup::setMasses(const std::vector<float>& newMasses) {masses = newMasses;}

void ParticleGroup::addParticle(const Vector3& position, const Vector3& velocity, float mass)
{
    positions.push_back(position);
    velocities.push_back(velocity);
    accelerations.emplace_back();
    masses.push_back(mass);
    nbParticles++;
}

float ParticleGroup::kineticEnergy() const
{
    float energy = 0.0f;
    for (int i = 0; i < nbParticles; i++)
    {
        const Vector3& velocity = velocities[i];
        float velocitySquared = velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z;
        energy += 0.5f * masses[i] * velocitySquared;
    }
    return energy;
}

void ParticleGroup::massCenter(Vector3& position, Vector3& velocity) const
{
    float totalMass = 0.0f;
    position = Vector3();
    velocity = Vector3();
    for (int i = 0; i < nbParticles; i++)
    {
        float mass = masses[i];
        totalMass += mass;
        position.x += positions[i].x * mass;
        position.y += positions[i].y * mass;
        position.z += positions[i].z * mass;
        velocity.x += velocities[i].x * mass;
        velocity.y += velocities[i].y * mass;
        velocity.z += velocities[i].z * mass;
    }
    if (totalMass > 0.0f)
    {
        position.x /= totalMass;
        position.y /= totalMass;
        position.z /= totalMass;
        velocity.x /= totalMass;
        velocity.y /= totalMass;
        velocity.z /= totalMass;
    }
}

ParticleGroup ParticleGroup::groupToCpu() const
{
    if (device == "cpu"){return *this;}
    ParticleGroup group(nbParticles, "cpu");
    group.positions = positions;
    group.velocities = velocities;
    group.accelerations = accelerations;
    group.masses = masses;
    return group;
}

ParticleGroup ParticleGroup::groupToGpu() const
{
    if (device == "gpu") {return *this;}
    ParticleGroup group(nbParticles, "gpu");
    group.positions = positions;
    group.velocities = velocities;
    group.accelerations = accelerations;
    group.masses = masses;
    return group;
}

int ParticleGroup::getNbParticles() const {return nbParticles;}
const std::vector<Vector3>& ParticleGroup::getPositions() const {return positions;}
const std::vector<Vector3>& ParticleGroup::getVelocities() const {return velocities;}
const std::vector<float>& ParticleGroup::getMasses() const {return masses;}
const std::string& ParticleGroup::getDevice() const {return device;}


void ParticleGroup::copyPositionsTo(float* out, int count) const
{
    const int n = std::min(count, nbParticles);
    for (int i = 0; i < n; ++i) {
        out[3*i+0] = positions[i].x;
        out[3*i+1] = positions[i].y;
        out[3*i+2] = positions[i].z;
    }
}

void ParticleGroup::copyVelocitiesTo(float* out, int count) const
{
    const int n = std::min(count, nbParticles);
    for (int i = 0; i < n; ++i) {
        out[3*i+0] = velocities[i].x;
        out[3*i+1] = velocities[i].y;
        out[3*i+2] = velocities[i].z;
    }
}

void ParticleGroup::copyMassesTo(float* out, int count) const
{
    const int n = std::min(count, nbParticles);
    for (int i = 0; i < n; ++i)
        out[i] = masses[i];
}