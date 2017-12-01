#pragma once

#include "interaction.hpp"



class AngularLennardJones
    : public Interaction
{
public:
    virtual void setup() override;
    virtual float translation(const Particle&, const Particle&) const override;
    virtual float rotation(const Particle&, const Particle&) const override;
    virtual float constrained(const Particle&, const Particle&) const override;
    virtual float constrainedOld(const Particle&, const Particle&) const override;
    virtual cartesian translation_force(const Particle&, const Particle&) const override;
    virtual cartesian rotation_force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
    
protected:
    float chi(const Particle&, const Particle&) const; 

    float kappa {0};
    float a {0};
    float b {0};
    float c {0};
};