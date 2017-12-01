#include "verlet.hpp"



void Verlet::step(const unsigned long& steps)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);
        // save forces
        for(auto& target : *target_range)
        {
            assert(target);
            target->save();
        }

        updateCoords();
        updateForces();
        updateVelocities();
    }
}



void Verlet::updateCoords()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    const auto dt = getParameters().dt;
    const auto dt2half = dt*dt*0.5f;

    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->setCoords( (target->coordsOld() + target->velocityOld()*dt + target->forceOld()*dt2half)/target->getMass());
    });
}



void Verlet::updateForces()
{
    vesDEBUG(__PRETTY_FUNCTION__)

    // first set to 0
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [](auto& target) 
    {
        assert(target);
        target->clearForce();
    });

    // second do calculation
    tbb::parallel_for(std::size_t(0),target_range->size(), std::size_t(1), [&](const std::size_t& i) 
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            const auto& target1 = target_range->operator[](i);
            const auto& target2 = target_range->operator[](j);
            assert(i!=j);
            assert(target1);
            assert(target2);
            assert(interaction);
            const Particle::cartesian translation_force = interaction->translation_force(target1,target2);
            target1->addForce(translation_force);
            target2->addForce((-1.f)*translation_force);
        }
    });
}



void Verlet::updateVelocities()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    const auto dt_half = 0.5f*getParameters().dt;
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->addVelocity( (target->forceOld() + target->force())*dt_half/target->getMass() );
    });
}