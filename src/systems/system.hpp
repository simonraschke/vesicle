#pragma once

#include <memory>
#include <vector>
#include <algorithm>
#include <tbb/cache_aligned_allocator.h>
#include "box.hpp"
#include "../particles/generator_interface.hpp"



class System
    : Box<PERIODIC::ON>
{
public:
    void clear();

    template<typename T>
    void addParticles(ParticleGenerator<T>&&);

protected:
    using Box<PERIODIC::ON>::distance;

private:
    // std::unique_ptr<AlgorithmInterface> algorithm {nullptr};
    std::vector<std::unique_ptr<ParticleInterface>, tbb::cache_aligned_allocator<ParticleInterface>> particles {};
};



template<typename T>
void System::addParticles(ParticleGenerator<T>&& gen)
{
    particles.reserve(particles.size()+gen.size());
    for(auto&& p : gen)
    {
        particles.emplace_back(p);
    }
    // std::move(gen.begin(),gen.end(), std::back_inserter(particles));
}