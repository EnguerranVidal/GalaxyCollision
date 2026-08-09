#include "bindings.h"

#include "calculators/Calculator.h"
#include "calculators/NewtonCalculator.h"
#include "calculators/BarnesHutCalculator.h"


void bindCalculator(py::module_& m)
{
    py::class_<Calculator>(m, "Calculator");
    py::class_<NewtonCalculator, Calculator>(m, "NewtonCalculator")
        .def(py::init<float>(), py::arg("gravitationalConstant") = 1.0f)
        .def("computeAccelerations", &NewtonCalculator::computeAccelerations,py::arg("particles"));

    py::class_<BarnesHutCalculator, Calculator>(m, "BarnesHutCalculator")
        .def(py::init<float, float>(), py::arg("theta") = 0.5f, py::arg("gravitationalConstant") = 1.0f)
        .def("computeAccelerations", &BarnesHutCalculator::computeAccelerations,py::arg("particles"));
}