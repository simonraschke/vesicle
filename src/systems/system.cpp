#include "system.hpp"



void System::clear()
{
    algorithm.reset(nullptr);
    particles.clear();
}



void System::startSimulation()
{
    int i = 0;
    while(true)
    {
        algorithm->step();
        std::cout << i++ << std::endl;
        std::cout << distance( *particles[0], *particles[1]) << std::endl;
        trajectory_writer->write();

        if(i>=1000) break;
    }
}



// void System::setParameters(Parameters prms)
// {
//     parameters = std::make_unique<Parameters>(prms);
// }