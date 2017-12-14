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

#include "andersen.hpp"



void AndersenThermostat::apply()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    // coordinates
    {
        auto& target = (*target_range)[enhance::random<std::size_t>(0,target_range->size()-1)];
        float mean = 0.f;
        float sigma = std::sqrt(getParameters().temperature/target->getMass());
        std::normal_distribution<float> gaussian(mean, sigma);
        float new_velocity = gaussian(enhance::RandomEngine.pseudo_engine);
        auto dimension = enhance::random<unsigned int>(0,2);
        auto velocities = target->velocity();
        velocities(dimension) = new_velocity;
        target->setVelocity(velocities);
    }

    // // orientation
    // {
    //     auto& target = (*target_range)[enhance::random<std::size_t>(0,target_range->size()-1)];
    //     float mean = 0.f;
    //     float sigma = std::sqrt(getParameters().temperature/target->getMass());
    //     std::normal_distribution<float> gaussian(mean, sigma);
    //     float new_velocity = gaussian(enhance::RandomEngine.pseudo_engine);
    //     auto dimension = enhance::random<unsigned int>(0,2);
    //     auto circularVelocities = target->circularVelocity();
    //     circularVelocities(dimension) = new_velocity;
    //     target->setCircularVelocity(circularVelocities);
    // }
}