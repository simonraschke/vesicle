#include "particle_base.hpp"



const ParticleInterface::cartesian& ParticleInterface::coords() const
{
    return *currentCoords;
}



const ParticleInterface::cartesian& ParticleInterface::orientation() const
{
    return *currentOrientation;
}



const ParticleInterface::cartesian& ParticleInterface::coordsOld() const
{
    return *oldCoords;
}



const ParticleInterface::cartesian& ParticleInterface::orientationOld() const
{
    return *oldOrientation;
}