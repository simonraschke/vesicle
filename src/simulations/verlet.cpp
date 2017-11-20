#include "verlet.hpp"



void Verlet::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);
        // std::cout << "target range" << std::endl;
        for(auto& target : *target_range)
        {
            assert(target);
            // std::cout << "target " << target->coords() << std::endl << std::endl;
            
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