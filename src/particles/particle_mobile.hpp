#pragma once

#include "particle_base.hpp"



class ParticleMobile 
    : public ParticleInterface
{
    virtual void updateCoords(cartesian&&) override;
    virtual void updateOrientation(cartesian&&) override;
};