#include "bindings.h"

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "particles/ParticleGroup.h"
#include "math/Vector3.h"

#include <stdexcept>
#include <algorithm>

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
        .def("getNbParticles", &ParticleGroup::getNbParticles)
        .def("getDevice", &ParticleGroup::getDevice, py::return_value_policy::reference_internal)
        .def("getPositions", &ParticleGroup::getPositions)
        .def("getVelocities", &ParticleGroup::getVelocities)
        .def("getMasses", &ParticleGroup::getMasses)
        .def("copyPositionsTo", [](const ParticleGroup& particles, py::array_t<float, py::array::c_style | py::array::forcecast> out)
            {
                py::buffer_info buf = out.request();
                const int n = particles.getNbParticles();
                if (buf.ndim == 2)
                {
                    if (buf.shape[0] < n || buf.shape[1] != 3)
                        throw std::runtime_error("copyPositionsTo: expected shape (N, 3)");
                }
                else if (buf.ndim == 1)
                {
                    if (buf.shape[0] < 3 * n)
                        throw std::runtime_error("copyPositionsTo: expected length >= 3*N");
                }
                else
                {
                    throw std::runtime_error("copyPositionsTo: array must be 1D or 2D");
                }

                particles.copyPositionsTo(static_cast<float*>(buf.ptr), n);
            },
            py::arg("out"))
        .def("copyVelocitiesTo", [](const ParticleGroup& particles, py::array_t<float, py::array::c_style | py::array::forcecast> out)
            {
                py::buffer_info buf = out.request();
                const int n = particles.getNbParticles();
                if (buf.ndim == 2)
                {
                    if (buf.shape[0] < n || buf.shape[1] != 3)
                        throw std::runtime_error("copyVelocitiesTo: expected shape (N, 3)");
                }
                else if (buf.ndim == 1)
                {
                    if (buf.shape[0] < 3 * n)
                        throw std::runtime_error("copyVelocitiesTo: expected length >= 3*N");
                }
                else
                {
                    throw std::runtime_error("copyVelocitiesTo: array must be 1D or 2D");
                }
                particles.copyVelocitiesTo(static_cast<float*>(buf.ptr), n);
            },
            py::arg("out"))
        .def("copyMassesTo", [](const ParticleGroup& particles, py::array_t<float, py::array::c_style | py::array::forcecast> out)
            {
                py::buffer_info buf = out.request();
                const int n = particles.getNbParticles();
                if (buf.ndim != 1 || buf.shape[0] < n)
                    throw std::runtime_error("copyMassesTo: expected 1D array of length >= N");
                particles.copyMassesTo(static_cast<float*>(buf.ptr), n);
            },
            py::arg("out"))
        .def("massCenter", [](const ParticleGroup& particles)
            {
                Vector3 position;
                Vector3 velocity;
                particles.massCenter(position, velocity);
                return py::make_tuple(position, velocity);
            });
}