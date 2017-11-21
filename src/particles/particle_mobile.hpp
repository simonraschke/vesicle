#pragma once

#include "particle.hpp"



class ParticleMobile 
    : public Particle
{
public:
    virtual void setCoords(const cartesian&) override;
    virtual void setVelocity(const cartesian&) override;
    virtual void setForce(const cartesian&) override;
    virtual void setOrientation(const cartesian&) override;

    virtual std::string name() const override;
};