#include "langevin.hpp"



void Langevin::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);
        for(auto& target : *target_range)
        {
            assert(target);
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



void Langevin::updateVelocities()
{

}


