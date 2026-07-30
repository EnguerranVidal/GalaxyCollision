#include "Calculator.h"


Calculator::Calculator(float gravitationalConstant): gravitationalConstant(gravitationalConstant), softening(1e-3f) {}

float Calculator::getGravitationalConstant() const {return gravitationalConstant;}

float Calculator::getSoftening() const {return softening;}