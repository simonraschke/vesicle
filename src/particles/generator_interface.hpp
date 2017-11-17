#pragma once

#include <vector>
#include "particle_mobile.hpp"
#include "particle_frame.hpp"



template<typename T, typename ENABLER = typename std::enable_if<std::is_base_of<ParticleInterface,T>::value>>
struct ParticleGenerator
{
    ParticleGenerator(std::size_t);

    std::size_t size() const;

    std::vector<ParticleInterface>::iterator begin() ;
    std::vector<ParticleInterface>::iterator end() ;

private:
    std::vector<std::unique_ptr<ParticleInterface>> particles;
};



template<typename T, typename ENABLER>
ParticleGenerator<T,ENABLER>::ParticleGenerator(std::size_t num)
    : particles(num)
{
    for(std::size_t i = 0; i < num; ++i)
    {
        particles.emplace_back(std::unique_ptr<ParticleInterface>( new T() ));
    }
}



template<typename T, typename ENABLER>
std::size_t ParticleGenerator<T,ENABLER>::size() const
{
    return particles.size();
}