#include "bindings.h"

#include "distributions/Distribution.h"
#include "distributions/BasicDistribution.h"
#include "distributions/GalaxyDistribution.h"


void bindDistribution(py::module_& m)
{
    py::class_<BasicDistribution, Distribution>(m, "BasicDistribution")
        .def(py::init<>())
        .def("generate", &BasicDistribution::generate, py::arg("parameters"));

    py::class_<GalaxyDistribution, Distribution>(m, "GalaxyDistribution")
        .def(py::init<>())
        .def("generate", &GalaxyDistribution::generate, py::arg("parameters"));
}