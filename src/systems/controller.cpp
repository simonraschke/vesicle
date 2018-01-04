/*  
*   Copyright 2017 Simon Raschke
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

#include "controller.hpp"



std::atomic<int> Controller::SIGNAL = {0};
tbb::mutex Controller::signal_mutex {};



void Controller::signal(int SIG)
{
    vesDEBUG(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
    tbb::mutex::scoped_lock lock(Controller::signal_mutex);
    SIGNAL.store(SIG);
}



void SimulationControl::make_nodes()
{
    start_node = std::make_unique<tbb::flow::broadcast_node<tbb::flow::continue_msg>>(flow);


    step_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            vesDEBUG("step_node")
            if(system.getAlgorithm())
                system.getAlgorithm()->step();
            else
                vesCRITICAL("no algorithm was set") 
        });


    thermostat_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            vesDEBUG("thermostat_node")
            if(system.getThermostat())
                system.getThermostat()->apply(); 
        });


    history_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            vesDEBUG("history_node")
            system.addTime(system.getParameters().dt); 
            HistoryBuffer buffer;
            buffer.time = std::make_unique<float>(system.getTime()); 
            // try
            // {
            //     tbb::parallel_invoke
            //     (
            //         [&]{ buffer.kineticEnergy = std::make_unique<float>(system.kineticEnergy()); },
            //         [&]{ buffer.potentialEnergy = std::make_unique<float>(system.potentialEnergy()); }
            //     );
            // }
            // catch(std::runtime_error& e){ vesWARNING(e.what()) }
            history_storage.flush(buffer);
        });


    trajectory_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            vesDEBUG("trajectory_node")
            if(system.getTrajectoryWriter())
                system.getTrajectoryWriter()->write(history_storage);
        });
}



void SimulationControl::setup()
{   
    vesDEBUG(__PRETTY_FUNCTION__)

    {
        // syytem is ParameterDependentComponent
        system.setParameters(getParameters());
    }

    // MUST
    {
        // add particles via factory class 
        system.addParticles(ParticleFactory<ParticleMobile>(getParameters().mobile));

        // and chose distribution
        system.distributeParticles<RandomDistributor>();
    }

    //MUST
    {
        // chose the algorithm depending on parameters
        if(getParameters().algorithm == std::string("verlet"))
            system.setAlgorithm<Verlet>();
        else if(getParameters().algorithm == "shakeVerlet" )
            system.setAlgorithm<ShakeVerlet>();
        else if(getParameters().algorithm == "langevin" )
            system.setAlgorithm<Langevin>();
        else if(getParameters().algorithm == "montecarlo" )
        {
            system.setAlgorithm<MonteCarlo>();

            // if algorithm is monte carlo add acceptance criterion
            if(getParameters().acceptance == std::string("metropolis"))
                system.getAlgorithm()->setAcceptance<Metropolis>();

            assert(system.getAlgorithm()->getAcceptance());
        }
        assert(system.getAlgorithm());
    }

    //MUST
    {
        // chose interactions
        if(getParameters().interaction == "lj")
            system.setInteraction<LennardJones>();
        else if(getParameters().interaction == "alj")
            system.setInteraction<AngularLennardJones>();
        assert(system.getInteraction());
    }

    // OPTIONAL
    {
        // add thermostat if wished
        // set nullptr if monte carlo because montecarlo implements temperature in acceptance
        if(getParameters().thermostat == "andersen")
        if(getParameters().algorithm  != "montecarlo")
        {
            system.setThermostat<AndersenThermostat>();
            assert(system.getThermostat());
        }
    }
    //OPTIONAL
    {   
        // add a trajectory writer if wished 
        if(getParameters().traj == "gro")
        {
            system.setTrajectoryWriter<TrajectoryWriterGro>();
            assert(system.getTrajectoryWriter());
        }

        // tell the trajectory writer if potential is anisotropic
        if(system.getTrajectoryWriter())
        {
            system.getTrajectoryWriter()->setAnisotropic(system.getInteraction()->isAnisotropic());
        }
    }
    
    // gernerate the node pointers
    make_nodes();

    // design of flow graph
    flow.reset();
    tbb::flow::make_edge(*start_node,*step_node);
    tbb::flow::make_edge(*step_node,*thermostat_node);
    tbb::flow::make_edge(*step_node,*history_node);
    tbb::flow::make_edge(*history_node,*trajectory_node);
}



void SimulationControl::start()
{
    vesDEBUG(__PRETTY_FUNCTION__)

    // print initial trajectory for start time
    {
        HistoryBuffer buffer;
        buffer.time = std::make_unique<float>(system.getTime()); 
        // try
        // {
        //     tbb::parallel_invoke
        //     (
        //         [&]{ buffer.kineticEnergy = std::make_unique<float>(system.kineticEnergy()); },
        //         [&]{ buffer.potentialEnergy = std::make_unique<float>(system.potentialEnergy()); }
        //     );
        // }
        // catch(std::runtime_error& e){ vesWARNING(e.what()) }
        history_storage.flush(buffer);
        if(system.getTrajectoryWriter())
            system.getTrajectoryWriter()->write(history_storage);
    }

    // step and time to track speed
    std::size_t i = 0;
    auto start = std::chrono::high_resolution_clock::now();

    while(SIGNAL.load() == 0)
    {
        // start flow graph from start_node
        start_node->try_put(tbb::flow::continue_msg());

        // and wait for it to finish all nodes
        flow.wait_for_all();

        // track speed
        if(i%system.getParameters().traj_skip==0) 
        {
            auto now = std::chrono::high_resolution_clock::now();
            auto dura = std::chrono::duration<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - start)).count();
            vesLOG("step " << std::setw(10) << std::right << i << " \t time:  " 
                << std::setw(10) << std::left << dura << " s \t time/step  " 
                << std::setw(10) << std::left << dura/system.getParameters().traj_skip*1000 << " ms \t time/step/particle  " 
                << std::setw(10) << std::left << dura/system.getParameters().traj_skip/system.getParticles().size()*1e6 << " ns")
            start = now;
        }
        ++i;

    }

    // write hstory on exit
    history_storage.dumpToFile("history.dat");
}



void SimulationControl::pause()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}