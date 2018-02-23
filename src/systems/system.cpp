/*  
*   Copyright 2017-2018 Simon Raschke
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

#include "system.hpp"



void System::clear()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    algorithm.reset(nullptr);
    thermostat.reset(nullptr);
    trajectory_writer.reset(nullptr);
    particles.clear();
    assert(!algorithm);
    assert(!thermostat);
    assert(!trajectory_writer);
}



const PARTICLERANGE& System::getParticles() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return particles;
}



PARTICLERANGE& System::getParticles()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return particles;
}



const std::unique_ptr<Algorithm>& System::getAlgorithm() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!algorithm) throw std::logic_error("algorithm is nullptr. set algorithm first");
    assert(algorithm);
    return algorithm;
}



const std::unique_ptr<Interaction>& System::getInteraction() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!algorithm) throw std::logic_error("algorithm is nullptr. set algorithm first");
    assert(algorithm);
    if(!algorithm->getInteraction()) throw std::logic_error("interaction is nullptr. set interaction first");
    assert(algorithm->getInteraction());
    return algorithm->getInteraction();
}



const std::unique_ptr<Thermostat>& System::getThermostat() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return thermostat;
}




const std::unique_ptr<TrajectoryWriter>& System::getTrajectoryWriter() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return trajectory_writer;
}



float System::potentialEnergy() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    std::atomic<float> energy_sum = {0.f};
    tbb::parallel_for(std::size_t(0),particles.size(), std::size_t(1), [&](const std::size_t& i)
    {
        float pre_sum = 0;
        for(std::size_t j = 0; j<i; ++j)
        {
            assert(algorithm);
            assert(algorithm->getInteraction());
            assert(particles[i]);
            assert(particles[j]);
            pre_sum += algorithm->getInteraction()->potential(particles[i],particles[j]) + algorithm->getInteraction()->potential(particles[i],particles[j]);
        }

        // this works for an atomic float, but may busy wait to compare_exchange like hell if too many particles
        auto current = energy_sum.load();
        while (!energy_sum.compare_exchange_weak(current, current + pre_sum))
            current = energy_sum.load();
        
    });
    return !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
}


float System::kineticEnergy() const
{   
    vesDEBUG(__PRETTY_FUNCTION__)
    float value = PARALLEL_REDUCE(float,particles,[&](float f, auto& p){ assert(p); return f + 0.5f*p->velocity().squaredNorm(); });
    return !std::isnan(value) ? value : throw std::runtime_error("kinetic Energy is NAN");
}



void System::addTime(float t)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    time_elapsed += t;
}



void System::setTime(float t)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    time_elapsed = t;
}



float System::getTime() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return time_elapsed;
}