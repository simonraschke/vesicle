#include "particle.hpp"



const Particle::cartesian& Particle::coords() const
{
    return *currentCoords;
}



const Particle::cartesian& Particle::orientation() const
{
    return *currentOrientation;
}



const Particle::cartesian& Particle::coordsOld() const
{
    return *oldCoords;
}



const Particle::cartesian& Particle::orientationOld() const
{
    return *oldOrientation;
}