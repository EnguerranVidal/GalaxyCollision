#include "bindings.h"

#include "integrators/EulerExplicitIntegrator.h"
#include "integrators/RK4Integrator.h"


void bindIntegrator(py::module_& m)
{
    py::class_<Integrator>(m, "Integrator");
    py::class_<EulerExplicitIntegrator>(m, "EulerIntegrator")
        .def(py::init<float>(), py::arg("timeStep"))
        .def("step", &EulerExplicitIntegrator::step, py::arg("particles"), py::arg("calculator"));

    py::class_<RK4Integrator>(m, "RK4Integrator")
        .def(py::init<float>(), py::arg("timeStep"))
        .def("step", &RK4Integrator::step, py::arg("particles"), py::arg("calculator"));
}