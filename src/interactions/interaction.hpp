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

    float translation(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    float rotation(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    cartesian translation_force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;
    cartesian rotation_force(const std::unique_ptr<Particle>&, const std::unique_ptr<Particle>&) const;

    virtual float translation(const Particle&, const Particle&) const = 0 ;
    virtual float rotation(const Particle&, const Particle&) const = 0 ;
    virtual cartesian translation_force(const Particle&, const Particle&) const = 0 ;
    virtual cartesian rotation_force(const Particle&, const Particle&) const = 0 ;


    virtual void setup();

    virtual bool isAnisotropic() const = 0;

protected:
    Interaction() = default;

private:
};



inline float Interaction::translation(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return translation(*ptr1,*ptr2);
}



inline float Interaction::rotation(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return rotation(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::translation_force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return translation_force(*ptr1,*ptr2);
}



inline Interaction::cartesian Interaction::rotation_force(const std::unique_ptr<Particle>& ptr1, const std::unique_ptr<Particle>& ptr2) const
{
    return rotation_force(*ptr1,*ptr2);
}



inline void Interaction::setup()
{
    return;
}