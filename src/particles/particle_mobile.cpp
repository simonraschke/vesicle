#include "particle_mobile.hpp"



void ParticleMobile::updateCoords(cartesian&& newCoords)
{
    assert(oldCoords);
    assert(currentCoords);
    std::swap(currentCoords,oldCoords);
    *currentCoords = std::move(newCoords);
}



void ParticleMobile::updateOrientation(cartesian&& newOrientation)
{
    assert(oldOrientation);
    assert(currentOrientation);
    std::swap(currentOrientation,oldOrientation);
    *currentOrientation = std::move(newOrientation);
}