#include <cstdlib>
#include <iostream>
#include "particles/particle_mobile.hpp"
#include "systems/system.hpp"


int main()
{
    System sys;
    sys.setParameters(Parameters());
    sys.addParticles(ParticleFactory<ParticleMobile>(2));
    sys.distributeParticles<RandomDistributor>();
    sys.setAlgorithm<Verlet>();
    sys.setInteraction<LennardJones>();
    sys.setTrajectoryWriter<TrajectoryWriterGro>();
    
    sys.startSimulation();

    return EXIT_SUCCESS;
}
