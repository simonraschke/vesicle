#pragma once

#include "particle.hpp"



class ParticleMobile 
    : public Particle
{
public:
    virtual void updateCoords(const cartesian&) override;
    virtual void updateOrientation(const cartesian&) override;
};