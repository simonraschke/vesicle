#include <cstdlib>
#include <iostream>
#include "particles/particle_mobile.hpp"
#include "systems/system.hpp"


int main()
{
    
    ParticleMobile p1;
    ParticleMobile p2;
    p1.updateCoords(Eigen::Vector3f(0,0,0));
    p2.updateCoords(Eigen::Vector3f(4,0,0));

    Box<PERIODIC::ON> box;
    box.setLengthX(10);
    box.setLengthY(10);
    box.setLengthZ(10);

    std::cout << box.distance(p1.coords(),p2.coords()) << std::endl;

    System sys;
    sys.addParticles(ParticleFactory<ParticleMobile>(5));

    return EXIT_SUCCESS;
}