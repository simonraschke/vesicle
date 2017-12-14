#pragma once

#include "interaction.hpp"



class LennardJones
    : public Interaction
{
    virtual float potential(const Particle&, const Particle&) const override;
    virtual cartesian force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
};