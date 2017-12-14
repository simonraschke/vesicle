/*  
*   Copyright 2017 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

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
            const Particle::cartesian translation_force = interaction->force(target1,target2);
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
    
}