#include "bindings.h"

#include "parameters/SolverParameters.h"
#include "parameters/BasicDistributionParameters.h"
#include "parameters/GalaxyDistributionParameters.h"


void bindParameters(py::module_& m)
{

    py::class_<BasicDistributionParameters>(m, "BasicDistributionParameters")
        .def(py::init<>())
        .def_readwrite("positionScale", &BasicDistributionParameters::positionScale)
        .def_readwrite("velocityScale", &BasicDistributionParameters::velocityScale)
        .def_readwrite("massMinimum", &BasicDistributionParameters::massMinimum)
        .def_readwrite("massMaximum", &BasicDistributionParameters::massMaximum);

    py::class_<GalaxyDistributionParameters>(m, "GalaxyDistributionParameters")
        .def(py::init<>())
        .def_readwrite("radius", &GalaxyDistributionParameters::radius)
        .def_readwrite("height", &GalaxyDistributionParameters::height)
        .def_readwrite("totalMass", &GalaxyDistributionParameters::totalMass)
        .def_readwrite("velocityDispersion", &GalaxyDistributionParameters::velocityDispersion)
        .def_readwrite("bulgeFraction", &GalaxyDistributionParameters::bulgeFraction)
        .def_readwrite("haloFraction", &GalaxyDistributionParameters::haloFraction)
        .def_readwrite("plummerRadius", &GalaxyDistributionParameters::plummerRadius)
        .def_readwrite("haloRadius", &GalaxyDistributionParameters::haloRadius);

    py::class_<SolverParameters>(m, "SolverParameters")
        .def(py::init<>())
        .def_readwrite("name", &SolverParameters::name)
        .def_readwrite("timeStep", &SolverParameters::timeStep)
        .def_readwrite("nbParticles", &SolverParameters::nbParticles)
        .def_readwrite("theta", &SolverParameters::theta)
        .def_readwrite("tileSize", &SolverParameters::tileSize)
        .def_readwrite("blockSize", &SolverParameters::blockSize)
        .def_readwrite("seed", &SolverParameters::seed)
        .def_readwrite("device", &SolverParameters::device)
        .def_readwrite("gravitationalConstant", &SolverParameters::gravitationalConstant)
        .def_readwrite("integratorType", &SolverParameters::integratorType)
        .def_readwrite("calculatorType", &SolverParameters::calculatorType)
        .def_readwrite("distributionType", &SolverParameters::distributionType)
        .def_readwrite("endless", &SolverParameters::endless)
        .def_readwrite("maxTime", &SolverParameters::maxTime)
        .def_readwrite("saveResults", &SolverParameters::saveResults)
        .def_readwrite("basicDistributionParameters", &SolverParameters::basicDistributionParameters)
        .def_readwrite("galaxyDistributionParameters", &SolverParameters::galaxyDistributionParameters);

}