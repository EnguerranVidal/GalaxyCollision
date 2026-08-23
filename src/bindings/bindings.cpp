#include "bindings.h"


PYBIND11_MODULE(engine, m)
{
    bindMath(m);
    bindParticles(m);
    bindParameters(m);
    bindCalculator(m);
    bindIntegrator(m);
}