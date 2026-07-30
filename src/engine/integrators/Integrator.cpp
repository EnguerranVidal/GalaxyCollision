#include "Integrator.h"


Integrator::Integrator(float timeStep): timeStep(timeStep) {}

float Integrator::getTimeStep() const {return timeStep;}