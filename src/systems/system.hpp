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

#pragma once

#include <memory>
#include <vector>
#include <algorithm>
// #include <tbb/cache_aligned_allocator.h>
#include "definitions.hpp"
#include "box.hpp"
#include "particles/particle_factory.hpp"
#include "particles/particle_distributor.hpp"
#include "simulations/verlet.hpp"
#include "simulations/shake_verlet.hpp"
#include "simulations/langevin.hpp"
#include "simulations/monte_carlo.hpp"
#include "interactions/lennard_jones.hpp"
#include "interactions/angular_lennard_jones.hpp"
#include "vesicleIO/parameters.hpp"
#include "vesicleIO/trajectory.hpp"
#include "thermostats/andersen.hpp"



class System
    : public Box<PERIODIC::ON>
    , public virtual ParameterDependentComponent
{
public:

    // reset the system
    void clear();

    // control
    template<typename T>
    void addParticles(ParticleFactory<T>&&);
    const PARTICLERANGE& getParticles() const;

    template<typename D,typename ENABLER = typename std::enable_if<std::is_base_of<Distributor,D>::value>::type>
    void distributeParticles();

    template<typename A,typename ENABLER = typename std::enable_if<std::is_base_of<Algorithm,A>::value>::type>
    void setAlgorithm();
    const std::unique_ptr<Algorithm>& getAlgorithm() const;

    template<typename I,typename ENABLER = typename std::enable_if<std::is_base_of<Interaction,I>::value>::type>
    void setInteraction();
    const std::unique_ptr<Interaction>& getInteraction() const;

    template<typename T,typename ENABLER = typename std::enable_if<std::is_base_of<Thermostat,T>::value>::type>
    void setThermostat();
    const std::unique_ptr<Thermostat>& getThermostat() const;

    template<typename W,typename ENABLER = typename std::enable_if<std::is_base_of<TrajectoryWriter,W>::value>::type>
    void setTrajectoryWriter();
    const std::unique_ptr<TrajectoryWriter>& getTrajectoryWriter() const;

    float kineticEnergy() const;
    float potentialEnergy() const;

    void addTime(float);
    float getTime() const;

    using ParameterDependentComponent::getParameters;
    
protected:
    using Box<PERIODIC::ON>::distance;
    using Box<PERIODIC::ON>::squared_distance;

private:
    PARTICLERANGE particles {};
    std::unique_ptr<Algorithm> algorithm {nullptr};
    std::unique_ptr<Thermostat> thermostat {nullptr};
    std::unique_ptr<TrajectoryWriter> trajectory_writer {nullptr};
    float time_elapsed {0.0};
};



template<typename T>
void System::addParticles(ParticleFactory<T>&& factory)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    particles.reserve(particles.size()+factory.size());
    while(factory)
    {
        particles.push_back(factory.createParticle());
    }
}



template<typename D,typename ENABLER>
void System::distributeParticles()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    D dist;
    dist.setParameters(getParameters());
    dist(&particles);
}



template<typename A,typename ENABLER>
void System::setAlgorithm()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    algorithm.reset(nullptr);
    assert(!algorithm);
    algorithm = std::make_unique<A>();
    assert(algorithm);
    algorithm->setParameters(getParameters());
    algorithm->setTarget(&particles);
    algorithm->setup();
}



template<typename T,typename ENABLER>
void System::setThermostat()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    thermostat.reset(nullptr);
    assert(!thermostat);
    thermostat = std::make_unique<T>();
    assert(thermostat);
    thermostat->setParameters(getParameters());
    thermostat->setTarget(&particles);
}



template<typename I,typename ENABLER>
void System::setInteraction()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(algorithm);
    algorithm->setInteraction<I>();
}



template<typename W,typename ENABLER>
void System::setTrajectoryWriter()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    trajectory_writer = std::make_unique<W>();
    assert(trajectory_writer);
    trajectory_writer->setParameters(getParameters());
    trajectory_writer->setFilename("trajectory");
    trajectory_writer->setTarget(&particles);
}