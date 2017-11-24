#include "particle.hpp"



void Particle::save()
{
    *oldCoords = *currentCoords;
    *oldVelocity = *currentVelocity;
    *oldForce = *currentForce;
    *oldOrientation = *currentOrientation;
    *oldCircularVelocity = *currentCircularVelocity;
    *oldTorque = *currentTorque;
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



void Particle::clearTorque()
{
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
    *currentVelocity += c;
}



void Particle::addForce(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentForce += c;
}



void Particle::addCircularVelocity(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentCircularVelocity += c;
}



void Particle::addTorque(const cartesian& c)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    *currentTorque += c;
}



// void Particle::addOrientation(const cartesian& c)
// {
//     tbb::spin_mutex::scoped_lock lock(mutex);
//     *currentOrientation += c;
// }



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



const Particle::cartesian& Particle::circularVelocity() const
{
    return *currentCircularVelocity;
}



const Particle::cartesian& Particle::circularVelocityOld() const
{
    return *oldCircularVelocity;
}



const Particle::cartesian& Particle::torque() const
{
    return *currentTorque;
}



const Particle::cartesian& Particle::torqueOld() const
{
    return *oldTorque;
}



void Particle::setMass(float m)
{
    mass = m;
}



float Particle::getMass() const
{
    return mass;
}