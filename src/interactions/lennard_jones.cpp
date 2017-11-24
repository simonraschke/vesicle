#include "lennard_jones.hpp"



float LennardJones::isotropic(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    return 4.f*(r6*r6-r6);
}



float LennardJones::anisotropic(const Particle& p1 __attribute__((unused)) , const Particle& p2  __attribute__((unused))) const 
{
    return 0.f;
}



LennardJones::cartesian LennardJones::isotropic_force(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    const float value = -24.f*r2*r6*(r6*2-1.f);

    return distance_vector(p1,p2)*value;
}



LennardJones::cartesian LennardJones::anisotropic_force(const Particle& p1 __attribute__((unused)) , const Particle& p2  __attribute__((unused))) const 
{
    return cartesian::Zero();
}



bool LennardJones::isAnisotropic() const
{
    return false;
}