#pragma once

#include "interaction_interface.hpp"



class LennardJones
    : public Interaction
{
    virtual float value(const Particle&, const Particle&) const override;
    virtual float derivative(const Particle&, const Particle&) const override;

private:
    float power6_term_reciprocal(const Particle&, const Particle&) const;
    float power6_term(const Particle&, const Particle&) const;
};