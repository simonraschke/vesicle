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
    using Box<PERIODIC::ON>::squared_distance;

    virtual ~Interaction() = default;

    float value(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    virtual float value(const Particle&, const Particle&) const = 0 ;

    float derivative(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    virtual float derivative(const Particle&, const Particle&) const = 0 ;

protected:
    Interaction() = default;

private:
};



inline float Interaction::value(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return value(*ptr1,*ptr2);
}



inline float Interaction::derivative(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return derivative(*ptr1,*ptr2);
}