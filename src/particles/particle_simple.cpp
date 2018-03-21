/*  
*   Copyright 2017-2018 Simon Raschke
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

#include "particle_simple.hpp"



ParticleSimple::ParticleSimple(PARTICLETYPE t)
    : position(cartesian::Zero())
    , orientation(cartesian::Zero())
    , velocity(cartesian::Zero())
    , force(cartesian::Zero())
    , mass(0)
    , type(t)
    , ID(0)
    , parent(nullptr)
{

}



ParticleSimple::ParticleSimple(Particle* other)
    : position(other->coords())
    , orientation(other->orientation())
    , velocity(other->velocity())
    , force(other->force())
    , mass(other->getMass())
    , type(other->getType())
    , ID(other->ID)
    , parent(other)
{

}



ParticleSimple::ParticleSimple(ParticleSimple* other)
    : position(other->position)
    , orientation(other->orientation)
    , velocity(other->velocity)
    , force(other->force)
    , mass(other->mass)
    , type(other->type)
    , ID(other->ID)
    , parent(nullptr)
{

}



bool ParticleSimple::operator==(Particle* other) const
{
    assert(other);
    assert(parent);
    return parent ? parent.get() == other : throw std::runtime_error("cannot compare ParticleSimple without parent to Particle*");
}



bool ParticleSimple::operator==(const Particle* other) const
{
    assert(other);
    assert(parent);
    return parent ? parent.get() == other : throw std::runtime_error("cannot compare ParticleSimple without parent to Particle*");
}



bool ParticleSimple::operator==(const ParticleSimple& other) const
{
    return this == std::addressof(other);
}