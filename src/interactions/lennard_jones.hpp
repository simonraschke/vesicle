#pragma once

#include "interaction.hpp"



class LennardJones
    : public Interaction
{
    virtual float value(const Particle&, const Particle&) const override;
    virtual cartesian force(const Particle&, const Particle&) const override;

private:
    // float power6_term_reciprocal(const Particle&, const Particle&) const;
    // float power6_term(const Particle&, const Particle&) const;
};