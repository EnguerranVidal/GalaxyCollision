#pragma once

#include <pybind11/pybind11.h>


namespace py = pybind11;

void bindMath(py::module_& m);
void bindParticles(py::module_& m);
void bindParameters(py::module_& m);
void bindCalculator(py::module_& m);
void bindIntegrator(py::module_& m);
void bindDistribution(py::module_& m);