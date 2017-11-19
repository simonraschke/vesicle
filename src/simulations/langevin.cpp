#include "langevin.hpp"



void Langevin::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        for(auto& target : target_range)
        {
            std::cout << target->coords() << std::endl << std::endl;
        }
    }
}



void Langevin::updateForces()
{

}



void Langevin::updateCoords()
{

}



void Langevin::updateOrientations()
{

}


