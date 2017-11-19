#pragma once

#include "vesicleIO/parameters.hpp"
#include "particles/particle_base.hpp"
#include "systems/box.hpp"

#include <cmath>



struct Interaction
    : public ParameterDependentComponent
    , public virtual Box<PERIODIC::ON>
{
    // typedef ParticleInterface ParticleInterface;
    // typedef ParticleInterface::cartesian cartesian;

    using Box<PERIODIC::ON>::distance;

    virtual ~Interaction() = default;

    float value(const std::unique_ptr<ParticleInterface>&, const std::unique_ptr<ParticleInterface>&) const;
    virtual float value(const ParticleInterface&, const ParticleInterface&) const = 0 ;

    float derivative(const std::unique_ptr<ParticleInterface>&, const std::unique_ptr<ParticleInterface>&) const;
    virtual float derivative(const ParticleInterface&, const ParticleInterface&) const = 0 ;

protected:
    Interaction() = default;

private:
    // float     
};



inline float Interaction::value(const std::unique_ptr<ParticleInterface>& ptr1, const std::unique_ptr<ParticleInterface>& ptr2) const
{
    return value(*ptr1,*ptr2);
}



inline float Interaction::derivative(const std::unique_ptr<ParticleInterface>& ptr1, const std::unique_ptr<ParticleInterface>& ptr2) const
{
    return derivative(*ptr1,*ptr2);
}