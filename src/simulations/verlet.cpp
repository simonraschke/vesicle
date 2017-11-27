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
        // updateOrientations();
        updateForces();
        updateVelocities();
    }
}



void Verlet::updateCoords()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    const auto dt = getParameters().dt;
    const auto dt2half = dt*dt*0.5f;

    // std::for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->setCoords( (target->coordsOld() + target->velocityOld()*dt + target->forceOld()*dt2half)/target->getMass());
        target->setOrientation( (target->orientationOld() + target->circularVelocityOld()*dt + target->torqueOld()*dt2half)/target->getMass() );
        
        const auto torque_vector = target->orientation().cross(target->orientationOld());
        const float angle = enhance::absolute_angle(target->orientationOld(),target->orientation());
        
        // rotate circular velocity and force
        // according to change of orientation
        {
            // const float decay = std::exp(-getParameters().dt);
            const float decay = 1;
            Eigen::AngleAxisf circular_velocity_rotation( angle, torque_vector.normalized() );
            target->setCircularVelocity( (circular_velocity_rotation * target->circularVelocityOld())*decay );
            target->setTorque( (circular_velocity_rotation * target->torqueOld())*decay*decay );
            vesDEBUG( "angle " << angle << "  torque" << torque_vector.normalized().format(ROWFORMAT) << "  circular velocity "<< target->circularVelocity().format(ROWFORMAT))
        }
        // else
        // {
            // vesDEBUG("WARINING: NO TORQUE CALCULATED, torque: " << torque_vector.normalized().format(ROWFORMAT) << ", angle: " << angle)
        // }
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

            if( std::abs(isotropic_force.normalized().dot(anisotropic_force.normalized()))-1.f > 1e-3f) 
                throw std::runtime_error("direction of forces should be the same, but scalar product is: "+std::to_string(isotropic_force.normalized().dot(anisotropic_force.normalized())));

            const Particle::cartesian chi_force_cartesian = interaction->chi_force(target1,target2);
            vesDEBUG("chi_force " << chi_force_cartesian.format(ROWFORMAT))
            target1->addTorque(chi_force_cartesian);
            target2->addTorque((-1.f)*chi_force_cartesian);
        }
    });
}



void Verlet::updateVelocities()
{
    vesDEBUG(__PRETTY_FUNCTION__<<"text")
    const auto dt_half = 0.5f*getParameters().dt;
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->addVelocity( (target->forceOld() + target->force())*dt_half/target->getMass() );
        target->addCircularVelocity( (target->torqueOld() + target->torque())*dt_half/target->getMass() );
    });
}