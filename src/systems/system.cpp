#include "system.hpp"



void System::clear()
{
    algorithm.reset(nullptr);
    particles.clear();
}



void System::startSimulation()
{
    algorithm->step();
}



// void System::setParameters(Parameters prms)
// {
//     parameters = std::make_unique<Parameters>(prms);
// }