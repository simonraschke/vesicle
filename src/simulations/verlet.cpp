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

        updateCoords();
        // updateOrientations();
        updateForces();
        updateVelocities();
    }
}



void Verlet::updateCoords()
{
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        // std::cout << target->velocity().norm() << std::endl;
        target->setCoords( (target->coordsOld() + target->velocityOld()*getParameters().dt + target->force()*0.5f*getParameters().dt*getParameters().dt)/target->getMass());
    });
}



void Verlet::updateOrientations()
{
    std::for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    // tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        // std::cout << target->velocity().norm() << std::endl;
        // target->setOrientation( (target->orientation() + target->velocity()*getParameters().dt + target->force()*0.5f*getParameters().dt*getParameters().dt)/target->getMass());
        auto torque = target->orientationOld().cross(target->forceOld());
        // std::cout << "torque: " << torque.format(ROWFORMAT) << std::endl;
        const Eigen::AngleAxisf rotate(torque.norm()*getParameters().dt/target->getMass(), torque);
        target->setOrientation( rotate * target->orientationOld() );
    });
}



void Verlet::updateForces()
{
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
            assert(target1);
            assert(target2);
            assert(interaction);
            const Particle::cartesian force_cartesian = interaction->force(target1,target2);
            target1->addForce(force_cartesian);
            // std::cout << force_cartesian.norm() << std::endl;
            target2->addForce((-1.f)*force_cartesian);
        }
    });
}



void Verlet::updateVelocities()
{
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->addVelocity( (target->forceOld() + target->force())*0.5f*getParameters().dt/target->getMass() );
    });
}