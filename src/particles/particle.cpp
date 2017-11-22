#include "particle.hpp"



void Particle::save()
{
    *oldCoords = *currentCoords;
    *oldVelocity = *currentVelocity;
    *oldForce = *currentForce;
    *oldOrientation = *currentOrientation;
}



void Particle::clearCoords()
{
    *currentCoords = cartesian::Zero();
}



void Particle::clearVelocity()
{
    *currentVelocity = cartesian::Zero();
}



void Particle::clearForce()
{
    *currentForce = cartesian::Zero();
}



void Particle::clearOrientation()
{
    *currentOrientation = cartesian::Zero();
}



void Particle::addCoords(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentCoords += c;
}



void Particle::addVelocity(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentVelocity += c;
}



void Particle::addForce(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentForce += c;
}



void Particle::addOrientation(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentOrientation += c;
}



const Particle::cartesian& Particle::coords() const
{
    return *currentCoords;
}



const Particle::cartesian& Particle::coordsOld() const
{
    return *oldCoords;
}



const Particle::cartesian& Particle::velocity() const
{
    return *currentVelocity;
}



const Particle::cartesian& Particle::velocityOld() const
{
    return *oldVelocity;
}



const Particle::cartesian& Particle::force() const
{
    return *currentForce;
}



const Particle::cartesian& Particle::forceOld() const
{
    return *oldForce;
}



const Particle::cartesian& Particle::orientation() const
{
    return *currentOrientation;
}



const Particle::cartesian& Particle::orientationOld() const
{
    return *oldOrientation;
}