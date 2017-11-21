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
#include "simulations/langevin.hpp"
#include "interactions/lennard_jones.hpp"
#include "vesicleIO/parameters.hpp"
#include "vesicleIO/trajectory.hpp"



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

    template<typename D,typename ENABLER = typename std::enable_if<std::is_base_of<Distributor,D>::value>::type>
    void distributeParticles();

    template<typename A,typename ENABLER = typename std::enable_if<std::is_base_of<Algorithm,A>::value>::type>
    void setAlgorithm();

    template<typename I,typename ENABLER = typename std::enable_if<std::is_base_of<Interaction,I>::value>::type>
    void setInteraction();

    template<typename W,typename ENABLER = typename std::enable_if<std::is_base_of<TrajectoryWriter,W>::value>::type>
    void setTrajectoryWriter();
    TrajectoryWriter& getTrajectoryWriter() const;

    void startSimulation();

protected:
    using Box<PERIODIC::ON>::distance;
    using Box<PERIODIC::ON>::squared_distance;
    using ParameterDependentComponent::getParameters;

private:
    std::unique_ptr<Algorithm> algorithm {nullptr};
    std::unique_ptr<TrajectoryWriter> trajectory_writer {nullptr};
    PARTICLERANGE particles {};

};



template<typename T>
void System::addParticles(ParticleFactory<T>&& factory)
{
    particles.reserve(particles.size()+factory.size());
    while(factory)
    {
        particles.push_back(factory.createParticle());
    }
}



template<typename D,typename ENABLER>
void System::distributeParticles()
{
    D dist;
    dist.setParameters(getParameters());
    dist(&particles);
}



template<typename A,typename ENABLER>
void System::setAlgorithm()
{
    algorithm.reset(nullptr);
    assert(!algorithm);
    algorithm = std::make_unique<A>();
    assert(algorithm);
    algorithm->setParameters(getParameters());
    algorithm->setTarget(&particles);
}   



template<typename I,typename ENABLER>
void System::setInteraction()
{
    assert(algorithm);
    algorithm->setInteraction<I>();
}   



template<typename W,typename ENABLER>
void System::setTrajectoryWriter()
{
    trajectory_writer = std::make_unique<W>();
    assert(trajectory_writer);
    trajectory_writer->setParameters(getParameters());
    trajectory_writer->setFilename("trajectory");
    trajectory_writer->setTarget(&particles);
}   