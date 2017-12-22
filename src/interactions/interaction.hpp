/*  
*   Copyright 2017 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

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
    using Box<PERIODIC::ON>::distanceVector;
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