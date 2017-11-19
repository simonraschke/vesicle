#pragma once

#include "interaction_interface.hpp"



class LennardJones
    : public Interaction
{
    virtual float value(const ParticleInterface&, const ParticleInterface&) const override;
    virtual float derivative(const ParticleInterface&, const ParticleInterface&) const override;

private:
    float power6_term(const ParticleInterface&, const ParticleInterface&) const;
};