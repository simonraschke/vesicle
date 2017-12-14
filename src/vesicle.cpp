#include <cstdlib>
#include <iostream>
#include <csignal>
#include "vesicleIO/program_options.hpp"
#include "systems/controller.hpp"



int main(int argc, const char *argv[])
{
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


    SimulationControl control;
    {
        Parameters prms;
        prms.programOptions.read(argc,argv);
        prms.setup();
        control.setParameters(prms);
    }
    control.setup();

    Eigen::initParallel();
#ifndef NDEBUG
    tbb::task_arena limited(1);
    limited.execute([&]
    {
        control.start();
    });
#else
    control.start();
#endif

    return EXIT_SUCCESS;
}
