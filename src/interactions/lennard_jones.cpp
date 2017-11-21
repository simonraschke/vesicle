#include "lennard_jones.hpp"



float LennardJones::value(const Particle& p1, const Particle& p2) const 
{
    const float power6_reciprocal = power6_term_reciprocal(p1,p2);
    return 4.f*(power6_reciprocal*power6_reciprocal-power6_reciprocal);
}



LennardJones::cartesian LennardJones::force(const Particle& p1, const Particle& p2) const 
{
    const float power6 = power6_term(p1,p2);
    const float force = 24.f*(power6-2)/(power6*power6*distance(p1.coordsOld(),p2.coordsOld()));

    return distance_vector(p1.coordsOld(),p2.coordsOld())*force;
}



float LennardJones::power6_term(const Particle& p1, const Particle& p2) const
{
    return std::pow(squared_distance(p1.coordsOld(),p2.coordsOld()),3);
}



float  LennardJones::power6_term_reciprocal(const Particle& p1, const Particle& p2) const
{
    return std::pow(1.f/squared_distance(p1.coordsOld(),p2.coordsOld()),3);
}