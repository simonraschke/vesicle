#pragma once

#include <vector>
#include "particle_mobile.hpp"
#include "particle_frame.hpp"
#include "enhance/random.hpp"


// Can also be used as a temporary object
template<typename T, typename ENABLER = typename std::enable_if<std::is_base_of<ParticleInterface,T>::value>>
struct ParticleFactory
{
    ParticleFactory(std::size_t);

    std::size_t size() const;

    std::unique_ptr<ParticleInterface> createParticle();

    std::vector<std::unique_ptr<ParticleInterface>>::iterator begin() ;
    std::vector<std::unique_ptr<ParticleInterface>>::iterator end() ;

    std::vector<std::unique_ptr<ParticleInterface>> particles;
};



template<typename T, typename ENABLER>
ParticleFactory<T,ENABLER>::ParticleFactory(std::size_t num)
    : particles(num)
{
    particles.reserve(num);
    for(std::size_t i = 0; i < num; ++i)
    {
        particles.emplace_back(createParticle());
    }
}



template<typename T, typename ENABLER>
std::unique_ptr<ParticleInterface> ParticleFactory<T,ENABLER>::createParticle()
{
    return std::make_unique<T>();
}



template<typename T, typename ENABLER>
std::size_t ParticleFactory<T,ENABLER>::size() const
{
    return particles.size();
}



template<typename T, typename ENABLER>
std::vector<std::unique_ptr<ParticleInterface>>::iterator ParticleFactory<T,ENABLER>::begin()
{
    return particles.begin();
}



template<typename T, typename ENABLER>
std::vector<std::unique_ptr<ParticleInterface>>::iterator ParticleFactory<T,ENABLER>::end()
{
    return particles.end();
}