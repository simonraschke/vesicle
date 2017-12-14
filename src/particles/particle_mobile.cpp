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



std::string ParticleMobile::name() const
{
    return "MOBIL";
}