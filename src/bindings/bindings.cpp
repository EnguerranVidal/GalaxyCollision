#include "bindings.h"


PYBIND11_MODULE(engine, m)
{
    m.doc() = "C++ N-body simulation engine";
    bindMath(m);
    bindParticles(m);
    bindParameters(m);
    bindCalculator(m);
    bindIntegrator(m);
    bindDistribution(m);
}