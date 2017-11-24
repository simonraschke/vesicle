#pragma once

#include "vesicleIO/parameters.hpp"
#include "particles/particle.hpp"
#include "systems/box.hpp"

#include <cmath>



struct Interaction
    : public Box<PERIODIC::ON>
    , virtual public ParameterDependentComponent
{
    using Box<PERIODIC::ON>::distance;
    using Box<PERIODIC::ON>::distance_vector;
    using Box<PERIODIC::ON>::squared_distance;
    typedef Particle::cartesian cartesian;

    virtual ~Interaction() = default;

    float isotropic(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    float anisotropic(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    cartesian isotropic_force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    cartesian anisotropic_force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    cartesian chi_force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;

    virtual float isotropic(const Particle&, const Particle&) const = 0 ;
    virtual float anisotropic(const Particle&, const Particle&) const = 0 ;
    virtual cartesian isotropic_force(const Particle&, const Particle&) const = 0 ;
    virtual cartesian anisotropic_force(const Particle&, const Particle&) const = 0 ;
    virtual cartesian chi_force(const Particle&, const Particle&) const = 0 ;


    virtual void setup();

    virtual bool isAnisotropic() const = 0;

protected:
    Interaction() = default;

private:
};



inline float Interaction::isotropic(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return isotropic(*ptr1,*ptr2);
}



inline float Interaction::anisotropic(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return anisotropic(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::isotropic_force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return isotropic_force(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::anisotropic_force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return anisotropic_force(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::chi_force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return chi_force(*ptr1,*ptr2);
}



inline void Interaction::setup()
{
    return;
}