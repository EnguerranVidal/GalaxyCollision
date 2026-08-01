#pragma once


struct BasicDistributionParameters
{
    int nbParticles = 1000;
    float positionScale = 10.0f;
    float velocityScale = 1.0f;
    float massMinimum = 0.1f;
    float massMaximum = 1.0f;
};