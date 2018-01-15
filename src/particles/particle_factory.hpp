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

#pragma once

#include <vector>
#include "particle_mobile.hpp"
#include "particle_frame.hpp"
#include "enhance/random.hpp"


// Can also be used as a temporary object
// this class is move assignable
//
// can be called as often as max size was set and creats particles of given type
// class will be false if not all particles were created else true
template<typename T, typename ENABLER = typename std::enable_if<std::is_base_of<Particle,T>::value>>
struct ParticleFactory
{
    ParticleFactory(std::size_t);

    // false if particles to be created
    explicit operator bool() const;

    // will create std::unique_ptr<> of paticle of given type
    // will throw if called too often
    // check operator bool() before calling
    PARTICLERANGE::value_type createParticle();

    // get size information
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