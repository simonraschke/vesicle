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

#include "particle.hpp"



bool Particle::operator==(const Particle& other)
{
    return std::addressof(*this) == std::addressof(other);
}



void Particle::save()
{
    assert(oldCoords);
    assert(currentCoords);
    assert(oldVelocity);
    assert(currentVelocity);
    assert(oldForce);
    assert(currentForce);
    assert(oldOrientation);
    assert(currentOrientation);

    *oldCoords = *currentCoords;
    *oldVelocity = *currentVelocity;
    *oldForce = *currentForce;
    *oldOrientation = *currentOrientation;
}



void Particle::clearCoords()
{
    assert(currentCoords);
    *currentCoords = cartesian::Zero();
}



void Particle::clearVelocity()
{
    assert(currentVelocity);
    *currentVelocity = cartesian::Zero();
}



void Particle::clearForce()
{
    assert(currentForce);
    *currentForce = cartesian::Zero();
}



void Particle::clearOrientation()
{
    assert(currentOrientation);
    *currentOrientation = cartesian::Zero();
}



// void Particle::addCoords(const cartesian& c)
// {
//     tbb::spin_mutex::scoped_lock lock(mutex);
//     *currentCoords += c;
// }



void Particle::addVelocity(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    assert(currentVelocity);
    assert(!c.hasNaN());
    assert(!currentVelocity->hasNaN());
    *currentVelocity += c;
}



void Particle::addForce(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    if(c.hasNaN())
        vesCRITICAL(c.format(ROWFORMAT))
    assert(currentForce);
    assert(!c.hasNaN());
    assert(!currentForce->hasNaN());
    *currentForce += c;
}



// void Particle::addOrientation(const cartesian& c)
// {
//     tbb::spin_mutex::scoped_lock lock(mutex);
//     *currentOrientation += c;
// }



const Particle::cartesian& Particle::coords() const
{
    assert(currentCoords);
    assert(!currentCoords->hasNaN());
    return *currentCoords;
}



const Particle::cartesian& Particle::coordsOld() const
{
    assert(oldCoords);
    assert(!oldCoords->hasNaN());
    return *oldCoords;
}



const Particle::cartesian& Particle::velocity() const
{
    assert(currentVelocity);
    assert(!currentVelocity->hasNaN());
    return *currentVelocity;
}



const Particle::cartesian& Particle::velocityOld() const
{
    assert(oldVelocity);
    assert(!oldVelocity->hasNaN());
    return *oldVelocity;
}



const Particle::cartesian& Particle::force() const
{
    assert(currentForce);
    assert(!currentForce->hasNaN());
    return *currentForce;
}



const Particle::cartesian& Particle::forceOld() const
{
    assert(oldForce);
    assert(!oldForce->hasNaN());
    return *oldForce;
}



const Particle::cartesian& Particle::orientation() const
{
    assert(currentOrientation);
    assert(!currentOrientation->hasNaN());
    return *currentOrientation;
}



const Particle::cartesian& Particle::orientationOld() const
{
    assert(oldOrientation);
    assert(!oldOrientation->hasNaN());
    return *oldOrientation;
}



void Particle::setMass(float m)
{
    assert(!std::isnan(m));
    mass = m;
}



float Particle::getMass() const
{
    assert(!std::isnan(mass));
    return mass;
}