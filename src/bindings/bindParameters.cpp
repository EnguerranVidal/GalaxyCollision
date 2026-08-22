#include "bindings.h"

#include "parameters/SimulatorParameters.h"
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

    py::class_<SimulatorParameters>(m, "SimulatorParameters")
        .def(py::init<>())
        .def_readwrite("name", &SimulatorParameters::name)
        .def_readwrite("timeStep", &SimulatorParameters::timeStep)
        .def_readwrite("nbParticles", &SimulatorParameters::nbParticles)
        .def_readwrite("theta", &SimulatorParameters::theta)
        .def_readwrite("tileSize", &SimulatorParameters::tileSize)
        .def_readwrite("blockSize", &SimulatorParameters::blockSize)
        .def_readwrite("seed", &SimulatorParameters::seed)
        .def_readwrite("device", &SimulatorParameters::device)
        .def_readwrite("gravitationalConstant", &SimulatorParameters::gravitationalConstant)
        .def_readwrite("integratorType", &SimulatorParameters::integratorType)
        .def_readwrite("calculatorType", &SimulatorParameters::calculatorType)
        .def_readwrite("distributionType", &SimulatorParameters::distributionType)
        .def_readwrite("endless", &SimulatorParameters::endless)
        .def_readwrite("maxTime", &SimulatorParameters::maxTime)
        .def_readwrite("saveResults", &SimulatorParameters::saveResults)
        .def_readwrite("basicDistributionParameters", &SimulatorParameters::basicDistributionParameters)
        .def_readwrite("galaxyDistributionParameters", &SimulatorParameters::galaxyDistributionParameters);

}