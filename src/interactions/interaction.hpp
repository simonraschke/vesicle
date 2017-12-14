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

    float potential(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    cartesian force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;

    virtual float potential(const Particle&, const Particle&) const = 0 ;
    virtual cartesian force(const Particle&, const Particle&) const = 0 ;


    virtual void setup();

    virtual bool isAnisotropic() const = 0;

protected:
    Interaction() = default;

private:
};



inline float Interaction::potential(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return potential(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return force(*ptr1,*ptr2);
}



inline void Interaction::setup()
{
    return;
}