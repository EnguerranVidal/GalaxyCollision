#include "bindings.h"
#include <pybind11/stl.h>

#include "particles/ParticleGroup.h"
#include "math/Vector3.h"


void bindParticles(py::module_& m)
{
    py::class_<ParticleGroup>(m, "ParticleGroup")
        .def(py::init<int, const std::string&>(), py::arg("nbParticles"), py::arg("device") = "CPU")
        .def("setPositions", &ParticleGroup::setPositions, py::arg("positions"))
        .def("setVelocities", &ParticleGroup::setVelocities, py::arg("velocities"))
        .def("setMasses", &ParticleGroup::setMasses, py::arg("masses"))
        .def("addParticle", &ParticleGroup::addParticle, py::arg("position"), py::arg("velocity"), py::arg("mass"))
        .def("kineticEnergy", &ParticleGroup::kineticEnergy)
        .def("groupToCpu", &ParticleGroup::groupToCpu)
        .def("groupToGpu", &ParticleGroup::groupToGpu)
        .def("getNbParticles",  &ParticleGroup::getNbParticles)
        .def("getDevice", &ParticleGroup::getDevice, py::return_value_policy::reference_internal)
        .def("getPositions", [](const ParticleGroup& p) { return p.positions; })
        .def("getVelocities", [](const ParticleGroup& p) { return p.velocities; })
        .def("getMasses", [](const ParticleGroup& p) { return p.masses; })
        .def("massCenter", [](const ParticleGroup& particles)
        {
            Vector3 position;
            Vector3 velocity;
            particles.massCenter(position, velocity);
            return py::make_tuple(position, velocity);
        }
        );


}