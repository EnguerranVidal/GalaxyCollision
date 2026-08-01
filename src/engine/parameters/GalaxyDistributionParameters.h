#pragma once


struct GalaxyDistributionParameters
{
    int nbParticles = 5000;
    float totalMass = 1000.0f;
    float radius = 15.0f;
    float height = 2.0f;
    float bulgeFraction = 0.2f;
    float haloFraction = 0.2f;
    float velocityDispersion = 0.1f;
    float plummerRadius = 3.0f;
    float haloRadius = 5.0f;
};