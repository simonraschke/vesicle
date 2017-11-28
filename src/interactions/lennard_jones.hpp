#pragma once

#include "interaction.hpp"



class LennardJones
    : public Interaction
{
    virtual float translation(const Particle&, const Particle&) const override;
    virtual float rotation(const Particle&, const Particle&) const override;
    virtual cartesian translation_force(const Particle&, const Particle&) const override;
    virtual cartesian rotation_force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
};