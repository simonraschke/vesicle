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

#include <vector>
#include "particle_mobile.hpp"
#include "particle_frame.hpp"
#include "enhance/random.hpp"


// Can also be used as a temporary object
template<typename T, typename ENABLER = typename std::enable_if<std::is_base_of<Particle,T>::value>>
struct ParticleFactory
{
    ParticleFactory(std::size_t);

    explicit operator bool() const;
    PARTICLERANGE::value_type createParticle();
    std::size_t size() const;

private:
    std::size_t num_to_create;
};



template<typename T, typename ENABLER>
ParticleFactory<T,ENABLER>::ParticleFactory(std::size_t num)
    : num_to_create(num)
{

}



template<typename T, typename ENABLER>
PARTICLERANGE::value_type ParticleFactory<T,ENABLER>::createParticle()
{
    if(num_to_create)
    {
        --num_to_create;
        return std::make_unique<T>();
    }
    else 
    {
        throw std::logic_error("cannot create more Particle instances. Factory empty");
    }
}



template<typename T, typename ENABLER>
std::size_t ParticleFactory<T,ENABLER>::size() const
{
    return num_to_create;
}



template<typename T, typename ENABLER>
ParticleFactory<T,ENABLER>::operator bool() const
{
    return num_to_create > 0;
}