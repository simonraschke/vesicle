#pragma once

#include "interaction.hpp"



class LennardJones
    : public Interaction
{
    virtual float isotropic(const Particle&, const Particle&) const override;
    virtual float anisotropic(const Particle&, const Particle&) const override;
    virtual cartesian isotropic_force(const Particle&, const Particle&) const override;
    virtual cartesian anisotropic_force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
};