#include "controller.hpp"



std::atomic<int> Controller::SIGNAL = {0};



void Controller::signal(int SIG)
{
    SIGNAL.store(SIG);
}



void SimulationControl::setup()
{   
    flow.reset();

    system.setParameters(Parameters());
    system.addParticles(ParticleFactory<ParticleMobile>(2));
    system.distributeParticles<RandomDistributor>();
    system.setAlgorithm<Verlet>();
    system.setInteraction<LennardJones>();
    system.setTrajectoryWriter<TrajectoryWriterGro>();

    start_node = std::make_unique<tbb::flow::broadcast_node<tbb::flow::continue_msg>>(flow);

    step_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg){ system.getAlgorithm().step(); });

    trajectory_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg){ system.getTrajectoryWriter().write(); });

    tbb::flow::make_edge(*start_node,*step_node);
    tbb::flow::make_edge(*step_node,*trajectory_node);
}



void SimulationControl::start()
{
    int i = 0;
    while(SIGNAL.load() == 0)
    {
        start_node->try_put(tbb::flow::continue_msg());
        flow.wait_for_all();
        std::cout << i++ << std::endl;
        if(i>=100000) break;
    }
}



void SimulationControl::pause()
{

}