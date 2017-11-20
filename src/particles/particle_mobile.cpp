#include "particle_mobile.hpp"



void ParticleMobile::updateCoords(const cartesian& newCoords)
{
    assert(oldCoords);
    assert(currentCoords);
    std::swap(currentCoords,oldCoords);
    *currentCoords = newCoords;
}



void ParticleMobile::updateOrientation(const cartesian& newOrientation)
{
    assert(oldOrientation);
    assert(currentOrientation);
    std::swap(currentOrientation,oldOrientation);
    *currentOrientation = newOrientation;
}