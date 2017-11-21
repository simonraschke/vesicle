#include "system.hpp"



void System::clear()
{
    algorithm.reset(nullptr);
    particles.clear();
}



Algorithm& System::getAlgorithm() const
{
    return *algorithm;
}




TrajectoryWriter& System::getTrajectoryWriter() const
{
    return *trajectory_writer;
}



void System::startSimulation()
{
}



// void System::setParameters(Parameters prms)
// {
//     parameters = std::make_unique<Parameters>(prms);
// }