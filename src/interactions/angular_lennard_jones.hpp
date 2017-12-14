#pragma once

#include "interaction.hpp"



class AngularLennardJones
    : public Interaction
{
public:
    virtual void setup() override;
    
    virtual float potential(const Particle&, const Particle&) const override;
    virtual cartesian force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
    
protected:
    float chi(const Particle&, const Particle&) const; 

    float kappa {0};
    float a {0};
    float b {0};
    float c {0};
};