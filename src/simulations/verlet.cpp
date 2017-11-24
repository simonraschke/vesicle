#include "verlet.hpp"



void Verlet::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);
        // save forces
        for(auto& target : *target_range)
        {
            // std::cout << "coords  : " << target->coords().format(ROWFORMAT) << std::endl;
            // std::cout << "velocity: " << target->velocity().format(ROWFORMAT) << std::endl;
            // std::cout << "force   : " << target->force().format(ROWFORMAT) << std::endl;
            // std::cout << "orientat: " << target->orientation().format(ROWFORMAT) << std::endl;
            // std::cout << "cricvelo: " << target->circularVelocity().format(ROWFORMAT) << std::endl;
            // std::cout << "torque  : " << target->torque().format(ROWFORMAT) << std::endl << std::endl;
            assert(target);
            target->save();
        }


        // tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target1) 
        // {
        //     target1->anisotropic_scalar = std::acumulate(target_range->begin(), target_range->end(), [&]( float f, auto& target2) { return interaction->anisotropic_force(target1,target2).norm(); }) / target.range();
        // });

        updateCoords();
        // updateOrientations();
        updateForces();
        updateVelocities();
    }
}



void Verlet::updateCoords()
{
    const auto dt = getParameters().dt;
    const auto dt2half = dt*dt*0.5f;


    // std::cout << "angle111 = " << enhance::rad_to_deg(enhance::angle(Eigen::Vector3f::UnitY(),Eigen::Vector3f::UnitX())) << std::endl;

    std::for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    // tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->setCoords( (target->coordsOld() + target->velocityOld()*dt + target->forceOld()*dt2half)/target->getMass());
        target->setOrientation( (target->orientationOld() + target->circularVelocityOld()*dt + target->torqueOld()*dt2half)/target->getMass() );

        // std::cout << "orien: " << target->orientation().cross(target->orientationOld()).format(ROWFORMAT) << std::endl;
        // std::cout << "orien: " << target->orientation().format(ROWFORMAT) << std::endl;
        // std::cout << "angle = " << enhance::rad_to_deg(enhance::angle(target->orientation(),target->orientationOld())) << std::endl;

        const auto torque = target->orientation().cross(target->orientationOld());
        if( torque.norm() > 0.00001 )
        {
            Eigen::AngleAxisf rotation (enhance::angle(target->orientation(),target->orientationOld()), target->orientation().cross(target->orientationOld()));
            target->setCircularVelocity( rotation * target->circularVelocityOld() );

        }

        // const auto torque = target->orientationOld().cross(target->torqueOld()*dt2half/target->getMass());
        // const auto torque = ((target->orientationOld() + target->circularVelocityOld()*dt + target->torqueOld()*dt2half)/target->getMass()).cross(target->torqueOld());
        // const auto torque = target->orientationOld().cross(target->torqueOld());
        // std::cout << target->torqueOld().format(ROWFORMAT) << "    " << torque.norm() << std::endl;
        // const Eigen::AngleAxisf rotate( torque.norm(), torque.normalized());
        // std::cout << target->circularVelocity().norm() << "  " << torque.norm() << std::endl;
        // target->setOrientation( rotate * target->orientationOld() );
        // target->setCircularVelocity( rotate * target->circularVelocityOld() );
    });
}



void Verlet::updateForces()
{
    // first set to 0
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [](auto& target) 
    {
        assert(target);
        target->clearForce();
        target->clearTorque();
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
            const Particle::cartesian isotropic_force = interaction->isotropic_force(target1,target2);
            const Particle::cartesian anisotropic_force = interaction->anisotropic_force(target1,target2);
            target1->addForce(isotropic_force + anisotropic_force);
            target2->addForce((-1.f)*(isotropic_force + anisotropic_force));

            const float dot_prod = isotropic_force.normalized().dot(anisotropic_force.normalized());
            if( std::abs(dot_prod)-1.f > 0.001f) throw std::runtime_error("direction of forces should be the same, but scalar product is: "+std::to_string(dot_prod));

            const Particle::cartesian chi_force_cartesian = interaction->chi_force(target1,target2);
            target1->addTorque(chi_force_cartesian);
            target2->addTorque((-1.f)*chi_force_cartesian);
        }
    });
}



void Verlet::updateVelocities()
{
    const auto dt_half = 0.5f*getParameters().dt;
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->addVelocity( (target->forceOld() + target->force())*dt_half/target->getMass() );
        target->addCircularVelocity( (target->torqueOld() + target->torque())*dt_half/target->getMass() );
    });
}