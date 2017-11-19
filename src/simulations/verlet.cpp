#include "verlet.hpp"



void Verlet::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        for(auto& target : target_range)
        {
            std::cout << target->coords() << std::endl << std::endl;
        }
    }
}



void Verlet::updateForces()
{
    
}



void Verlet::updateCoords()
{

}



void Verlet::updateOrientations()
{

}