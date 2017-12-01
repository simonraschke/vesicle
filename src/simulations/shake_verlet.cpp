#include "shake_verlet.hpp"



void ShakeVerlet::step(const unsigned long& steps)
{
    vesDEBUG(__PRETTY_FUNCTION__)

    if(!jacobian) jacobian =           std::make_unique<Eigen::MatrixXf>(Eigen::MatrixXf::Zero(target_range->size(),target_range->size()));
    if(!jacobianOld) jacobianOld =     std::make_unique<Eigen::MatrixXf>(Eigen::MatrixXf::Zero(target_range->size(),target_range->size()));
    if(!lagrangian) lagrangian =       std::make_unique<Eigen::MatrixXf>(Eigen::MatrixXf::Zero(target_range->size(),target_range->size()));
    if(!lagrangianOld) lagrangianOld = std::make_unique<Eigen::MatrixXf>(Eigen::MatrixXf::Zero(target_range->size(),target_range->size()));

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



void ShakeVerlet::updateCoords()
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



void ShakeVerlet::updateForces()
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

    shake();
}



void ShakeVerlet::updateVelocities()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    const auto dt_half = 0.5f*getParameters().dt;
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        target->addVelocity( (target->forceOld() + target->force())*dt_half/target->getMass() );
    });
}



void ShakeVerlet::shake()
{
    const float optimum = 0.0;
    const float tolerance = 0.5;
    assert(target_range);

    jacobian.swap(jacobianOld);
    lagrangian.swap(lagrangianOld);

    // calculate all constraints
    // for(std::size_t i = 0; i < target_range->size(); ++i)
    tbb::parallel_for(std::size_t(0),target_range->size(), std::size_t(1), [&](const std::size_t& i) 
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            assert(i!=j);
            const auto& target1 = target_range->operator[](i);
            const auto& target2 = target_range->operator[](j);
            const float constrained = interaction->constrained(target1,target2);
            (*jacobian)(i,j) = constrained > 1e-3 ? constrained : 0.f;
            // jacobian(j,i) = constrained > 1e-3 ? constrained : 0.f;
        }
    });

    // if( jacobian->maxCoeff() <= tolerance )
        // return;
    
    tbb::parallel_for(std::size_t(0),target_range->size(), std::size_t(1), [&](const std::size_t& i) 
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            assert(i!=j);
            const auto& target1 = target_range->operator[](i);
            const auto& target2 = target_range->operator[](j);

            const float jacobianVal = (*jacobian)(i,j);
            const float jacobianOldVal = (*jacobianOld)(i,j);
            const float mass1 = 1.f/target1->getMass();
            const float mass2 = 1.f/target2->getMass();

            if( ((*jacobian)(i,j)) < 1e-15 ) continue;
            else if( jacobianOldVal < 1e-15 ) continue;
            else if( (jacobianVal*jacobianVal - optimum*optimum) < 1e-15 ) continue;

            const float lagrangianVal = (*lagrangian)(i,j);
            const float lagrangianOldVal = (*lagrangianOld)(i,j);

            (*lagrangian)(i,j) = (jacobianVal*jacobianVal - optimum*optimum) / 
                ( 2.f*(mass1 + mass2) * jacobianVal * jacobianOldVal );
            
            Eigen::Matrix<float,3,1> orientation_temp = target1->orientation();
            float add = lagrangianVal*lagrangianOldVal;
            orientation_temp = orientation_temp.unaryExpr([&](float f){ return f + add;});
            target1->setOrientation(orientation_temp); 
        }
    });

    auto constraints = (*jacobian)*(*lagrangian);

    // tbb::parallel_for(std::size_t(0),target_range->size(), std::size_t(1), [&](const std::size_t& i) 
    // {
    //     for(std::size_t j = 0; j<i; ++j)
    //     {
    //         assert(i!=j);
    //         const auto& target1 = target_range->operator[](i);
    //         const auto& target2 = target_range->operator[](j);

    //         const float constrainedOld = interaction->constrainedOld(*target1,*target2);


    //     }
    // });


    vesDEBUG('\n' << *jacobian)
    vesDEBUG('\n' << *lagrangian)
    vesDEBUG('\n' << constraints)
    // std::terminate();
}