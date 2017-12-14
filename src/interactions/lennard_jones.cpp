#include "lennard_jones.hpp"



float LennardJones::potential(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    return 4.f*(r6*r6-r6);
}



LennardJones::cartesian LennardJones::force(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1.coordsOld(),p2.coordsOld());
    const float r6 = r2*r2*r2;
    const float value = -24.f*r2*r6*(r6*2-1.f);
#ifndef NDEBUG
    if(!std::isfinite(squared_distance(p1,p2)) || !std::isfinite(r2) || !std::isfinite(value) || distance_vector(p1,p2).hasNaN()) 
    {
        vesWARNING("distance_vector " << distance_vector(p1,p2).format(ROWFORMAT))
        vesWARNING("squared_distance(p1,p2) " << squared_distance(p1,p2))
        vesWARNING("r2 " << r2)
        vesWARNING("r6 " << r6)
        vesWARNING("value " << value)
    }
#endif

    return distance_vector(p1.coordsOld(),p2.coordsOld())*value;
}



bool LennardJones::isAnisotropic() const
{
    return false;
}