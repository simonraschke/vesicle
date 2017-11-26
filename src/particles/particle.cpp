#include "particle.hpp"



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
    assert(oldCircularVelocity);
    assert(currentCircularVelocity);
    assert(oldTorque);
    assert(currentTorque);

    *oldCoords = *currentCoords;
    *oldVelocity = *currentVelocity;
    *oldForce = *currentForce;
    *oldOrientation = *currentOrientation;
    *oldCircularVelocity = *currentCircularVelocity;
    *oldTorque = *currentTorque;
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



void Particle::clearTorque()
{
    assert(currentTorque);
    *currentTorque = cartesian::Zero();
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
    assert(currentForce);
    assert(!c.hasNaN());
    assert(!currentForce->hasNaN());
    *currentForce += c;
}



void Particle::addCircularVelocity(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    assert(currentCircularVelocity);
    assert(!c.hasNaN());
    assert(!currentCircularVelocity->hasNaN());
    *currentCircularVelocity += c;
}



void Particle::addTorque(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    assert(currentTorque);
    assert(!c.hasNaN());
    assert(!currentTorque->hasNaN());
    *currentTorque += c;
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



const Particle::cartesian& Particle::circularVelocity() const
{
    assert(currentCircularVelocity);
    assert(!currentCircularVelocity->hasNaN());
    return *currentCircularVelocity;
}



const Particle::cartesian& Particle::circularVelocityOld() const
{
    assert(oldCircularVelocity);
    assert(!oldCircularVelocity->hasNaN());
    return *oldCircularVelocity;
}



const Particle::cartesian& Particle::torque() const
{
    assert(currentTorque);
    assert(!currentTorque->hasNaN());
    return *currentTorque;
}



const Particle::cartesian& Particle::torqueOld() const
{
    assert(oldTorque);
    assert(!oldTorque->hasNaN());
    return *oldTorque;
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