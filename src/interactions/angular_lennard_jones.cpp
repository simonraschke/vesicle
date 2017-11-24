#include "angular_lennard_jones.hpp"


#include <iostream>
float AngularLennardJones::isotropic(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    return 4.f*(r6*r6-r6);
}



AngularLennardJones::cartesian AngularLennardJones::isotropic_force(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    const float r6 = r2*r2*r2;
    return distance_vector(p1,p2)*((-24.f)*r2*r6*(r6*2-1.f));
}



float AngularLennardJones::anisotropic(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    return 4.f*chi(p1,p2)*r2*r2*r2;
}



AngularLennardJones::cartesian AngularLennardJones::anisotropic_force(const Particle& p1, const Particle& p2) const 
{
    const float r2 = 1.f/squared_distance(p1,p2);
    float r8 = r2*r2;
    r8 *= r8;
    return distance_vector(p1,p2)*(24.f*chi(p1,p2)*r8);
}



AngularLennardJones::cartesian AngularLennardJones::chi_force(const Particle& p1, const Particle& p2) const 
{
    const cartesian p1_orien_kappa = p1.orientation()*kappa/2.f;
    const cartesian p2_orien_kappa = p2.orientation()*kappa/2.f;
    const cartesian partA = distance_vector(p1_orien_kappa,p2_orien_kappa)*(distance_vector(p1_orien_kappa,p2_orien_kappa).norm()-a)*2;
    const cartesian partB = distance_vector(-p1_orien_kappa,-p2_orien_kappa)*(distance_vector(-p1_orien_kappa,-p2_orien_kappa).norm()-b)*2;
    const cartesian partC1 = distance_vector(p1_orien_kappa,-p2_orien_kappa)*(distance_vector(p1_orien_kappa,-p2_orien_kappa).norm()-c)*2;
    const cartesian partC2 = distance_vector(-p1_orien_kappa,p2_orien_kappa)*(distance_vector(-p1_orien_kappa,p2_orien_kappa).norm()-c)*2;

    return (partA + partB + partC1 + partC2)*(-1.f);
}



void AngularLennardJones::setup()
{
    kappa = getParameters().kappa;
    a = 1.f + kappa*std::sin(getParameters().gamma);
    b = 1.f - kappa*std::sin(getParameters().gamma);
    c = (cartesian(b,0,0) + cartesian(a-1.f, kappa*std::cos(getParameters().gamma), 0)).norm();
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



bool AngularLennardJones::isAnisotropic() const
{
    return true;
}
// float LennardJones::power6_term(const Particle& p1, const Particle& p2) const
// {
//     return std::pow(squared_distance(p1.coordsOld(),p2.coordsOld()),3);
// }



// float  LennardJones::power6_term_reciprocal(const Particle& p1, const Particle& p2) const
// {
//     return std::pow(1.f/squared_distance(p1.coordsOld(),p2.coordsOld()),3);
// }