#include "controller.hpp"



std::atomic<int> Controller::SIGNAL = {0};
tbb::mutex Controller::signal_mutex {};



void Controller::signal(int SIG)
{
    tbb::mutex::scoped_lock lock(Controller::signal_mutex);
    SIGNAL.store(SIG);
}



void SimulationControl::setup()
{   
    flow.reset();

    system.setParameters(Parameters());
    system.addParticles(ParticleFactory<ParticleMobile>(2));
    system.distributeParticles<RandomDistributor>();
    system.setAlgorithm<Verlet>();
    system.setThermostat<AndersenThermostat>();
    system.setInteraction<AngularLennardJones>();
    system.setTrajectoryWriter<TrajectoryWriterGro>();
    system.getTrajectoryWriter().setSkip(system.getParameters().trajectory_skip);
    system.getTrajectoryWriter().setAnisotropic(system.getAlgorithm().getInteraction().isAnisotropic());

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
    // tbb::flow::make_edge(*step_node,*thermostat_node);
    tbb::flow::make_edge(*step_node,*history_node);
    tbb::flow::make_edge(*history_node,*trajectory_node);
    

    // // write first trajectory step
    // {
    //     HistoryBuffer buffer;
    //     buffer.time = std::make_unique<float>(system.getTime()); 
    //     try
    //     {
    //         tbb::parallel_invoke
    //         (
    //             [&]{ buffer.kineticEnergy = std::make_unique<float>(system.kineticEnergy()); },
    //             [&]{ buffer.potentialEnergy = std::make_unique<float>(system.potentialEnergy()); }
    //         );
    //     }
    //     catch(std::runtime_error e){ std::cout << e.what() << std::endl;  }
    //     history_storage.flush(buffer);
    //     system.getTrajectoryWriter().write(history_storage);
    // }
}



void SimulationControl::start()
{
    history_node->try_put(tbb::flow::continue_msg());
    flow.wait_for_all();

    int i = 0;
    while(SIGNAL.load() == 0)
    {
        start_node->try_put(tbb::flow::continue_msg());
        flow.wait_for_all();
        if(i%system.getParameters().trajectory_skip==0) std::cout << i << std::endl;
        if(i++>=10000) break;
    }
    history_storage.dumpToFile("history.dat");
}



void SimulationControl::pause()
{

}