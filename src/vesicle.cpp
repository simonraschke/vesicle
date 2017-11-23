#include <cstdlib>
#include <iostream>
#include <csignal>
#include "systems/controller.hpp"


int main()
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
    control.setup();
    control.start();

    // std::cout << sizeof(ParticleMobile)*8 << std::endl;
    // AngularLennardJones interaction;
    // interaction.setParameters(Parameters());
    // interaction.setup();
    // ParticleMobile p1;
    // ParticleMobile p2;

    // p1.setCoords(Eigen::Vector3f(1,1,1));
    // p2.setCoords(Eigen::Vector3f(1.9,1,1));
    // p1.setOrientation(Eigen::Vector3f(+1.0,0,0));
    // p2.setOrientation(Eigen::Vector3f(-1.0,0,0));
    // for( int i = 0; i < 500; ++i)
    // {
    //     p1.save();
    //     p2.save();
    //     p2.setCoords(p2.coordsOld()+Eigen::Vector3f(0.01,0,0));
    //     std::cout << 0.9+i*0.01 << "\t" <<  interaction.value(p1,p2) << "\t" << interaction.force(p1,p2).norm() << std::endl;
    // }

    return EXIT_SUCCESS;
}
