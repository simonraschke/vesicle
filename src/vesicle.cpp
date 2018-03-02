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

#include "vesicleIO/parameters.hpp"
#include "systems/controller.hpp"
#include <csignal>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>



int main(int argc, const char *argv[]) try
{
    // register important signals in Controller bas class
    // allowing civilized shutdown
    std::signal( SIGHUP,  Controller::signal );
    std::signal( SIGINT,  Controller::signal );
    std::signal( SIGQUIT, Controller::signal );
    std::signal( SIGILL,  Controller::signal );
    std::signal( SIGTRAP, Controller::signal );
    std::signal( SIGABRT, Controller::signal );
    std::signal( SIGIOT,  Controller::signal );
    std::signal( SIGBUS,  Controller::signal );
    std::signal( SIGFPE,  Controller::signal );
    std::signal( SIGKILL, Controller::signal );
    std::signal( SIGTERM, Controller::signal );

    Parameters prms;
    prms.read(argc,argv);
    prms.setup();

    // create a task arena 
#ifndef NDEBUG
    tbb::task_scheduler_init init(1);
    tbb::task_arena limited(1);
#else
    tbb::task_scheduler_init init();
    tbb::task_arena limited(prms.cpu_threads);
#endif
    
    // execute in limited task arena
    limited.execute([&]
    {
        SimulationControl control;
        control.setParameters(prms);
        control.setup();
        control.start();
    });

    return EXIT_SUCCESS;
}
catch(std::exception& e)
{
    vesLOG("caught exception at main scope: " << e.what() )
}
