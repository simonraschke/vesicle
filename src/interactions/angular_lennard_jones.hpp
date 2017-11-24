#pragma once

#include "interaction.hpp"



class AngularLennardJones
    : public Interaction
{
// public:
    virtual void setup() override;
    virtual float isotropic(const Particle&, const Particle&) const override;
    virtual float anisotropic(const Particle&, const Particle&) const override;
    virtual cartesian isotropic_force(const Particle&, const Particle&) const override;
    virtual cartesian anisotropic_force(const Particle&, const Particle&) const override;

    virtual cartesian chi_force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
    
protected:
    float chi(const Particle&, const Particle&) const; 

    float kappa {0};
    float a {0};
    float b {0};
    float c {0};
    // float power6_term_reciprocal(const Particle&, const Particle&) const;
    // float power6_term(const Particle&, const Particle&) const;
};