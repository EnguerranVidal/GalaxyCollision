#pragma once

#include <string>

#include "BasicDistributionParameters.h"
#include "GalaxyDistributionParameters.h"


struct SolverParameters
{
    std::string name = "Untitled Simulation";
    float timeStep = 0.001f;
    int nbParticles = 1000;
    float theta = 0.5f;
    int tileSize = 16;
    int blockSize = 256;
    int seed = -1;
    std::string device = "GPU";
    float gravitationalConstant = 1.0f;
    std::string integratorType = "RK4";
    std::string distributionType = "GALAXY";
    std::string calculatorType = "BARNES_HUT";
    bool endless = true;
    float maxTime = 1000.0f;
    bool saveResults = false;
    BasicDistributionParameters basicDistributionParameters;
    GalaxyDistributionParameters galaxyDistributionParameters;
};