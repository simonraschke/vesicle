#include <cstdlib>
#include <iostream>
#include "particles/particle_mobile.hpp"
#include "systems/system.hpp"


int main()
{
    System sys;
    sys.setParameters(Parameters());
    sys.addParticles(ParticleFactory<ParticleMobile>(100));
    sys.distributeParticles<RandomDistributor>();
    sys.setAlgorithm<Verlet>();
    
    sys.startSimulation();

    return EXIT_SUCCESS;
}