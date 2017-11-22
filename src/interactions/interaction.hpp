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

    float value(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    virtual float value(const Particle&, const Particle&) const = 0 ;

    cartesian force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    virtual cartesian force(const Particle&, const Particle&) const = 0 ;

protected:
    Interaction() = default;

private:
};



inline float Interaction::value(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return value(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return force(*ptr1,*ptr2);
}