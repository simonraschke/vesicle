#include "controller.hpp"



std::atomic<int> Controller::SIGNAL = {0};
tbb::mutex Controller::signal_mutex {};



void Controller::signal(int SIG)
{
    vesDEBUG(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
    tbb::mutex::scoped_lock lock(Controller::signal_mutex);
    SIGNAL.store(SIG);
}



void SimulationControl::setup()
{   
    vesDEBUG(__PRETTY_FUNCTION__)
    flow.reset();

    system.setParameters(Parameters());
    system.addParticles(ParticleFactory<ParticleMobile>(20));
    system.distributeParticles<RandomDistributor>();
    system.setAlgorithm<Verlet>();
    system.setThermostat<AndersenThermostat>();
    system.setInteraction<AngularLennardJones>();
    system.setTrajectoryWriter<TrajectoryWriterGro>();
    system.getTrajectoryWriter().setSkip(system.getParameters().trajectory_skip);
    system.getTrajectoryWriter().setAnisotropic(system.getAlgorithm().getInteraction().isAnisotropic());

    // system.getParticles()[0]->setCoords(Eigen::Vector3f(1,1,1));
    // system.getParticles()[1]->setCoords(Eigen::Vector3f(2.32246204831,1,1));
    // system.getParticles()[0]->setOrientation(Eigen::Vector3f(0,1,0));
    // system.getParticles()[1]->setOrientation(Eigen::Vector3f(0,1,0));

    start_node = std::make_unique<tbb::flow::broadcast_node<tbb::flow::continue_msg>>(flow);

    step_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg){ system.getAlgorithm().step(); });

    thermostat_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg){ system.getThermostat().apply(); });

    history_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            system.addTime(system.getParameters().dt); 
            HistoryBuffer buffer;
            buffer.time = std::make_unique<float>(system.getTime()); 
            try
            {
                tbb::parallel_invoke
                (
                    [&]{ buffer.kineticEnergy = std::make_unique<float>(system.kineticEnergy()); },
                    [&]{ buffer.potentialEnergy = std::make_unique<float>(system.potentialEnergy()); }
                );
            }
            catch(std::runtime_error e){ std::cout << e.what() << std::endl;  }
            history_storage.flush(buffer);
        });

    trajectory_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg){ system.getTrajectoryWriter().write(history_storage); });

    tbb::flow::make_edge(*start_node,*step_node);
    tbb::flow::make_edge(*step_node,*thermostat_node);
    tbb::flow::make_edge(*step_node,*history_node);
    tbb::flow::make_edge(*history_node,*trajectory_node);
}



void SimulationControl::start()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    history_node->try_put(tbb::flow::continue_msg());
    flow.wait_for_all();

    int i = 0;
    while(SIGNAL.load() == 0)
    {
        start_node->try_put(tbb::flow::continue_msg());
        flow.wait_for_all();
        if(i%system.getParameters().trajectory_skip==0) std::cout << i << std::endl;
        if(i++>=500000) break;
    }
    history_storage.dumpToFile("history.dat");

    // auto v1 = Eigen::Vector3f(-1,-1,0);
    // auto v2 = Eigen::Vector3f( 1,0,0);
    // vesDEBUG( enhance::rad_to_deg( enhance::absolute_angle(v1,v2) ) )
}



void SimulationControl::pause()
{
    vesDEBUG(__PRETTY_FUNCTION__)

}