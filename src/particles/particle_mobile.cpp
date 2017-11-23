#include "particle_mobile.hpp"



void ParticleMobile::setCoords(const cartesian& newCoords)
{
    assert(currentCoords);
    *currentCoords = newCoords;
}



void ParticleMobile::setVelocity(const cartesian& newVelocity)
{
    assert(currentVelocity);
    *currentVelocity = newVelocity;
}



void ParticleMobile::setForce(const cartesian& newForce)
{
    assert(currentForce);
    *currentForce = newForce;
}



void ParticleMobile::setOrientation(const cartesian& newOrientation)
{
    assert(currentOrientation);
    *currentOrientation = newOrientation.normalized();
}



void ParticleMobile::setCircularVelocity(const cartesian& newCircularVelocity)
{
    assert(currentCircularVelocity);
    *currentCircularVelocity = newCircularVelocity;
}



void ParticleMobile::setTorque(const cartesian& newTorque)
{
    assert(currentTorque);
    *currentTorque = newTorque;
}



std::string ParticleMobile::name() const
{
    return "MOBIL";
}