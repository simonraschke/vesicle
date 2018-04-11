/*  
*   Copyright 2017-2018 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#include "angular_lennard_jones.hpp"

float AngularLennardJones::potential(const Particle& p1, const Particle& p2) const 
{
    if(p1.getType() == OSMOTIC)
    {
        if(p2.getType() == OSMOTIC)
            return 0.f;
        else
            return osmoticPotential(p1,p2);
    }
    else if(p2.getType() == OSMOTIC)
    {
        if(p1.getType() == OSMOTIC)
            return 0.f;
        else
            return osmoticPotential(p1,p2);
    }
    else
    {
        cartesian distance_vec = distanceVector(p1, p2);

        const float r2 = sigma/distance_vec.squaredNorm();
        if(r2 < cutoff_rez_sq) return 0.f;

        const float r6 = r2*r2*r2;

        distance_vec.normalize();
        const cartesian p1_orien_kappa = p1.orientation()*kappa/2.f;
        const cartesian p2_orien_kappa = p2.orientation()*kappa/2.f;

        const float chi = 
            std::pow(cartesian( -p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - a,2)
            + std::pow(cartesian(  p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - b,2)
            + std::pow(cartesian( -p1_orien_kappa + distance_vec - p2_orien_kappa ).norm() - c,2)
            + std::pow(cartesian(  p1_orien_kappa + distance_vec + p2_orien_kappa ).norm() - c,2);
            
        return 4.f*epsilon*(r6*r6-(1.f-chi)*r6);
    }
}



float AngularLennardJones::osmoticPotential(const Particle& p1, const Particle& p2) const 
{
    assert(p1.getType() == OSMOTIC || p2.getType() == OSMOTIC);
    assert(p1.getType() != p2.getType());
    
    // const Particle& osmotic = p1.getType() == OSMOTIC ? p1 : p2;
    const Particle& nonosmotic = p1.getType() != OSMOTIC ? p1 : p2;
    assert(std::addressof(p1) != std::addressof(p2));

    // attractive part
    float attractive = 0;
    {
        const cartesian nonosmotic_orien_kappa = nonosmotic.orientation()*kappa/2.f;
        cartesian distance_vec = distanceVector(p1, p2) - nonosmotic_orien_kappa;

        const float r2 = sigma/distance_vec.squaredNorm();
        if(r2 < cutoff_rez_sq) return 0.f;

        const float r6 = r2*r2*r2;

        distance_vec.normalize();
        attractive = r6*r6-r6;
    }

    float repulsive = 0;
    {
        const cartesian nonosmotic_orien_kappa = nonosmotic.orientation()*kappa/2.f;
        cartesian distance_vec = distanceVector(p1, p2) + nonosmotic_orien_kappa;

        const float r2 = sigma/distance_vec.squaredNorm();
        if(r2 < cutoff_rez_sq) return 0.f;
        const float r6 = r2*r2*r2;
        repulsive = r6;
    }

    return attractive + 4.f*epsilon*repulsive;
}



AngularLennardJones::cartesian AngularLennardJones::force(const Particle& p1, const Particle& p2) const 
{
    const cartesian distance_vec = distanceVector(p1, p2);

    const float r2 = sigma/distance_vec.squaredNorm();
    const float r6 = r2*r2*r2;

    const auto normed_dist_vec = distance_vec.normalized();
    const cartesian p1_orien_kappa = p1.orientation()*kappa/2.f;
    const cartesian p2_orien_kappa = p2.orientation()*kappa/2.f;

    const float chi = 
           std::pow(cartesian( -p1_orien_kappa + normed_dist_vec + p2_orien_kappa ).norm() - a,2)
         + std::pow(cartesian(  p1_orien_kappa + normed_dist_vec - p2_orien_kappa ).norm() - b,2)
         + std::pow(cartesian( -p1_orien_kappa + normed_dist_vec - p2_orien_kappa ).norm() - c,2)
         + std::pow(cartesian(  p1_orien_kappa + normed_dist_vec + p2_orien_kappa ).norm() - c,2);

    return distance_vec*epsilon*((-24.f)*r2*r6*(r6*2 + chi - 1.f));
}



void AngularLennardJones::setup()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    kappa = getParameters().kappa;
    a = 1.f + kappa*std::sin(getParameters().gamma);
    b = 1.f - kappa*std::sin(getParameters().gamma);
    c = (cartesian(b,0,0) + cartesian(a-1.f, kappa*std::cos(getParameters().gamma), 0)).norm();
    epsilon = getParameters().LJepsilon;
    sigma = getParameters().LJsigma;
    cutoff_rez_sq = 1.f/getParameters().cell_min_edge;
    cutoff_rez_sq *= cutoff_rez_sq;
}




// float AngularLennardJones::chi(const Particle& p1, const Particle& p2) const 
// {
//     const cartesian normed_dist_vec = distanceVector(p1, p2).normalized() * 1.f;
//     const cartesian p1_orien_kappa = p1.orientation()*kappa/2.f;
//     const cartesian p2_orien_kappa = p2.orientation()*kappa/2.f;

//     return std::pow(cartesian( -p1_orien_kappa + normed_dist_vec + p2_orien_kappa ).norm() - a,2)
//          + std::pow(cartesian(  p1_orien_kappa + normed_dist_vec - p2_orien_kappa ).norm() - b,2)
//          + std::pow(cartesian( -p1_orien_kappa + normed_dist_vec - p2_orien_kappa ).norm() - c,2)
//          + std::pow(cartesian(  p1_orien_kappa + normed_dist_vec + p2_orien_kappa ).norm() - c,2);
// }



bool AngularLennardJones::isAnisotropic() const
{
    return true;
}