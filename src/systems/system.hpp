#pragma once

#include <memory>
#include <vector>
#include <algorithm>
#include <tbb/cache_aligned_allocator.h>
#include "box.hpp"
#include "particles/particle_factory.hpp"
#include "simulations/langevuin.hpp"



class System
    : Box<PERIODIC::ON>
{
public:
    void clear();

    template<typename T>
    void addParticles(ParticleFactory<T>&&);

    template<typename A>
    void setAlgorithm();

protected:
    using Box<PERIODIC::ON>::distance;

private:
    std::unique_ptr<Algorithm> algorithm {nullptr};
    // std::vector<std::unique_ptr<ParticleInterface>, tbb::cache_aligned_allocator<ParticleInterface>> particles {};
    std::vector<std::unique_ptr<ParticleInterface>> particles {};
};



template<typename T>
void System::addParticles(ParticleFactory<T>&& gen)
{
    particles.reserve(particles.size()+gen.size());
    std::move(gen.begin(),gen.end(), std::back_inserter(particles));
}