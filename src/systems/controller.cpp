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

#include "controller.hpp"



std::atomic<int> Controller::SIGNAL = {0};
tbb::mutex Controller::signal_mutex {};



void Controller::signal(int SIG)
{
    static std::size_t got_called = 0;
    ++got_called;
    
    if(got_called == 1)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
    }
    else if (got_called == 2)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
        vesWARNING(__PRETTY_FUNCTION__ << " still trying civilized shutdown...")
    }
    else if (got_called == 3)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
        vesWARNING(__PRETTY_FUNCTION__ << " TERMINATING!")
        std::exit(SIG);
    }

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


    system_adjust_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            vesDEBUG("system_adjust_node")
            system.addTime(system.getParameters().dt);
        });


    trajectory_node = std::make_unique<tbb::flow::continue_node<tbb::flow::continue_msg>>
        (flow, [&](tbb::flow::continue_msg)
        { 
            vesDEBUG("trajectory_node")
            if(system.getTrajectoryWriter())
                system.getTrajectoryWriter()->write(system.getTime());
        });
}



void SimulationControl::setup()
{   
    vesDEBUG(__PRETTY_FUNCTION__)
    GLOBAL::getInstance().status.store(GLOBAL::PREPARATION);

    {
        // syytem is ParameterDependentComponent
        system.setParameters(getParameters());
    }

    // MUST
    {
        // add particles via factory class 
        std::size_t guiding_elements_all = std::pow(getParameters().frame_guides_grid_edge, 3) * getParameters().guiding_elements_each;
        system.addParticles(ParticleFactory<ParticleFrame>(guiding_elements_all));
        system.addParticles(ParticleFactory<ParticleMobile>(getParameters().mobile));
        system.addParticles(ParticleFactory<ParticleOsmotic>(getParameters().osmotic));
        
        // for(const auto& p : system.getParticles())
        // {
        //     vesLOG(*p);
        // }
        vesLOG(getParameters().guiding_elements_each << " " <<  getParameters().mobile << " " << getParameters().osmotic);

        // and chose distribution
        if(GLOBAL::getInstance().mode == GLOBAL::NEWRUN)
        {
            if(getParameters().osmotic > 0)
            {
                system.distributeParticles<OsmoticSystemDistributor>();
            }
            else if (getParameters().frame_guides_grid_edge > 0 && getParameters().guiding_elements_each > 0)
            {
                system.distributeParticles<FrameGuidedGridDistributor>();
            }
            else if (getParameters().guiding_elements_plane > 0)
            {
                system.distributeParticles<FrameGuidedPlaneDistributor>();
            }
            else
            {
                system.distributeParticles<RandomDistributor>();
            }
        }
        else if(GLOBAL::getInstance().mode == GLOBAL::RESTART && getParameters().in_traj == std::string("gro"))
        {
            system.distributeParticles<TrajectoryDistributorGro>();
        }


        if(GLOBAL::getInstance().mode == GLOBAL::RESTART)
        {
            TrajectoryReaderGro reader;
            reader.setParameters(getParameters());
            reader.setPath(getParameters().in_traj_path);
            reader.readAllFrames();
            system.setTime(reader.getTime());
            vesLOG("set start time to " << system.getTime())
        }
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
        if(getParameters().out_traj == "gro")
        {
            system.setTrajectoryWriter<TrajectoryWriterGro>();
            assert(system.getTrajectoryWriter());
        }

        // tell the trajectory writer if potential is anisotropic
        assert(system.getTrajectoryWriter());
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
    tbb::flow::make_edge(*step_node,*system_adjust_node);
    tbb::flow::make_edge(*system_adjust_node,*trajectory_node);
}



void SimulationControl::start()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    GLOBAL::getInstance().status.store(GLOBAL::RUNNING);

    if(std::abs(system.getTime() - getParameters().time_max) < 1e-3 ) 
    {
        vesLOG("SYSTEM ALREADY REACHED MAX TIME")
        return;
    }

    // step and time to track speed
    std::size_t i = 1;
    double dura = 0;
    double since_beginning = 0;
    auto the_very_beginning = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::high_resolution_clock::now();

    // print initial trajectory for start time
    if(system.getTrajectoryWriter() && GLOBAL::getInstance().mode == GLOBAL::NEWRUN)
        system.getTrajectoryWriter()->write(system.getTime(),true);

    const auto time_max = getParameters().time_max;

    while(SIGNAL.load() == 0 && system.getTime() <= time_max)
    {
        // start flow graph from start_node
        start_node->try_put(tbb::flow::continue_msg());

        // and wait for it to finish all nodes
        flow.wait_for_all();

        // track speed
        if(i%system.getParameters().out_traj_skip==0) 
        {
            now = std::chrono::high_resolution_clock::now();
            dura = std::chrono::duration<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - start)).count();
            since_beginning = std::chrono::duration<double>(std::chrono::duration_cast<std::chrono::seconds>(now - the_very_beginning)).count();
            vesLOG( "simulation time elapsed " <<  std::setw(14) << std::right << std::setprecision(3) << std::fixed << system.getTime() << " | step " 
                << std::setw(10) << std::right << i << " | time:  " 
                << std::setw(10) << std::right << dura << " s | time/step  " 
                << std::setw(6) << std::right << dura/system.getParameters().out_traj_skip*1000 << " ms | time/step/particle  " 
                << std::setw(6) << std::right << dura/system.getParameters().out_traj_skip/system.getParticles().size()*1e6 << " ns | "
                << "time per day " << std::setw(12) << std::scientific << std::right << (24.0*60*60)/(since_beginning?since_beginning:1)*(getParameters().dt*i)
                << " | time to target " << std::setw(12) << std::fixed << std::right << getParameters().time_max/((24.0*60*60)/since_beginning*(getParameters().dt*i)) << " d"
            )
            start = now;
        }
        ++i;
    }
}



void SimulationControl::pause()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}