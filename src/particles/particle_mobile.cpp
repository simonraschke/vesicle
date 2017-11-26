#include "particle_mobile.hpp"



void ParticleMobile::setCoords(const cartesian& newCoords)
{
    assert(!newCoords.hasNaN());
    assert(currentCoords);
    *currentCoords = newCoords;
}



void ParticleMobile::setVelocity(const cartesian& newVelocity)
{
    assert(!newVelocity.hasNaN());
    assert(currentVelocity);
    *currentVelocity = newVelocity;
}



void ParticleMobile::setForce(const cartesian& newForce)
{
    assert(!newForce.hasNaN());
    assert(currentForce);
    *currentForce = newForce;
}



void ParticleMobile::setOrientation(const cartesian& newOrientation)
{
    assert(!newOrientation.hasNaN());
    assert(currentOrientation);
    *currentOrientation = newOrientation.normalized();
}



void ParticleMobile::setCircularVelocity(const cartesian& newCircularVelocity)
{
    assert(!newCircularVelocity.hasNaN());
    assert(currentCircularVelocity);
    *currentCircularVelocity = newCircularVelocity;
}



void ParticleMobile::setTorque(const cartesian& newTorque)
{
    assert(!newTorque.hasNaN());
    assert(currentTorque);
    *currentTorque = newTorque;
}



std::string ParticleMobile::name() const
{
    return "MOBIL";
}