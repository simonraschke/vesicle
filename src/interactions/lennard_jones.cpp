#include "lennard_jones.hpp"



float LennardJones::value(const ParticleInterface&, const ParticleInterface&) const 
{
    return 0;
}



float LennardJones::derivative(const ParticleInterface&, const ParticleInterface&) const 
{
    return 0;
}



float LennardJones::power6_term(const ParticleInterface& p1, const ParticleInterface& p2) const
{
    // return std::pow(1.f/this->distance(p1,p2),6);
    return 0;
}