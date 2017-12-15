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
            try
            {
                tbb::parallel_invoke
                (
                    [&]{ buffer.kineticEnergy = std::make_unique<float>(system.kineticEnergy()); },
                    [&]{ buffer.potentialEnergy = std::make_unique<float>(system.potentialEnergy()); }
                );
            }
            catch(std::runtime_error& e){ vesWARNING(e.what()) }
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
        system.setParameters(getParameters());
    }
    // MUST
    {
        system.addParticles(ParticleFactory<ParticleMobile>(getParameters().mobile));
        system.distributeParticles<RandomDistributor>();
    }
    //MUST
    {
        if(getParameters().algorithm == std::string("verlet"))
            system.setAlgorithm<Verlet>();
        else if(getParameters().algorithm == "shakeVerlet" )
            system.setAlgorithm<ShakeVerlet>();
        else if(getParameters().algorithm == "langevin" )
            system.setAlgorithm<Langevin>();
        else if(getParameters().algorithm == "montecarlo" )
        {
            system.setAlgorithm<MonteCarlo>();

            if(getParameters().acceptance == std::string("metropolis"))
                system.getAlgorithm()->setAcceptance<Metropolis>();

            assert(system.getAlgorithm()->getAcceptance());
        }
        assert(system.getAlgorithm());
    }
    //MUST
    {
        if(getParameters().interaction == "lj")
            system.setInteraction<LennardJones>();
        else if(getParameters().interaction == "alj")
            system.setInteraction<AngularLennardJones>();
        assert(system.getInteraction());
    }
    // OPTIONAL
    {
        if(getParameters().thermostat == "andersen")
        if(getParameters().algorithm  != "montecarlo")
        {
            system.setThermostat<AndersenThermostat>();
            assert(system.getThermostat());
        }
    }
    //OPTIONAL
    {
        if(getParameters().traj == "gro")
        {
            system.setTrajectoryWriter<TrajectoryWriterGro>();
            assert(system.getTrajectoryWriter());
        }
        if(system.getTrajectoryWriter())
        {
            system.getTrajectoryWriter()->setAnisotropic(system.getInteraction()->isAnisotropic());
        }
    }
    
    make_nodes();

    // make flow graph
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
        try
        {
            tbb::parallel_invoke
            (
                [&]{ buffer.kineticEnergy = std::make_unique<float>(system.kineticEnergy()); },
                [&]{ buffer.potentialEnergy = std::make_unique<float>(system.potentialEnergy()); }
            );
        }
        catch(std::runtime_error& e){ vesWARNING(e.what()) }
        history_storage.flush(buffer);
        if(system.getTrajectoryWriter())
            system.getTrajectoryWriter()->write(history_storage);
    }

    std::size_t i = 0;
    while(SIGNAL.load() == 0)
    {
        start_node->try_put(tbb::flow::continue_msg());
        flow.wait_for_all();
        if(i%system.getParameters().traj_skip==0) std::cout << i << std::endl;
        ++i;
    }
    history_storage.dumpToFile("history.dat");
}



void SimulationControl::pause()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}