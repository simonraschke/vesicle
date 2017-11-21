#include "verlet.hpp"



void Verlet::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);

        // save forces
        for(auto& target : *target_range)
        {
            assert(target);
            target->save();
        }

        // calculate new positions
        for(auto& target : *target_range)
        {
            assert(target);
            target->setCoords( target->coords() + target->velocity()*getParameters().dt + target->force()*0.5f*getParameters().dt*getParameters().dt);
        }

        // calculate new forces
        updateForces();

    }
}



void Verlet::updateForces()
{
    // first set to 0
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [](auto& target) 
    // for(auto& target : *target_range)
    {
        assert(target);
        target->clearForce();
    }
    );

    // second do calculation
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target1) 
    // for(auto& target1 : *target_range)
    {
        for(auto& target2 : *target_range)
        {
            if( &(*target1) == &(*target2)) continue;
            // std::cout << &(*target1) << "  " << &(*target2) << std::endl;
            assert(target1);
            assert(target2);
            assert(interaction);
            const Particle::cartesian force_cartesian = interaction->force(target1,target2);
            target1->addForce(force_cartesian);
            target2->addForce((-1.f)*force_cartesian);
        }
    }
    );
}



void Verlet::updateCoords()
{

}



void Verlet::updateOrientations()
{

}