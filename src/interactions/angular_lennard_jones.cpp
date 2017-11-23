#include "angular_lennard_jones.hpp"


#include <iostream>
float AngularLennardJones::value(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    // std::cout << "chi=" << chi(p1,p2) << std::endl;
    return 4.f*(r6*r6-(1.f-chi(p1,p2))*r6);
}



AngularLennardJones::cartesian AngularLennardJones::force(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    const float value = -24.f*r2*r6*((1.f-chi(p1,p2))-r6*2);

    return distance_vector(p1,p2)*value;
}



float AngularLennardJones::chi(const Particle& p1, const Particle& p2) const 
{
    const cartesian normed_dist_vec = distance_vector(p1, p2).normalized() * 1.f;
    const cartesian p1_orien_kappa = p1.orientation()*kappa/2.f;
    const cartesian p2_orien_kappa = p2.orientation()*kappa/2.f;
    return std::pow(cartesian( -p1_orien_kappa + normed_dist_vec + p2_orien_kappa ).norm() - a,2)
         + std::pow(cartesian(  p1_orien_kappa + normed_dist_vec - p2_orien_kappa ).norm() - b,2)
         + std::pow(cartesian( -p1_orien_kappa + normed_dist_vec - p2_orien_kappa ).norm() - c,2)
         + std::pow(cartesian(  p1_orien_kappa + normed_dist_vec + p2_orien_kappa ).norm() - c,2);
}


// #include <iostream>
void AngularLennardJones::setup()
{
    kappa = getParameters().kappa;
    a = 1.f + kappa*std::sin(getParameters().gamma);
    b = 1.f - kappa*std::sin(getParameters().gamma);
    c = (cartesian(b,0,0) + cartesian(a-1.f, kappa*std::cos(getParameters().gamma), 0)).norm();
}





// float LennardJones::power6_term(const Particle& p1, const Particle& p2) const
// {
//     return std::pow(squared_distance(p1.coordsOld(),p2.coordsOld()),3);
// }



// float  LennardJones::power6_term_reciprocal(const Particle& p1, const Particle& p2) const
// {
//     return std::pow(1.f/squared_distance(p1.coordsOld(),p2.coordsOld()),3);
// }