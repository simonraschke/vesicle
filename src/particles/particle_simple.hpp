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

#include "particle.hpp"
#include "enhance/observer_ptr.hpp"
#include "enhance/concurrent_container.hpp"
#include <tbb/atomic.h>
#include <eigen3/Eigen/Core>
#include <tbb/concurrent_vector.h>



class ParticleSimple
{
public:
    typedef Particle::cartesian cartesian;

    // construction
    explicit ParticleSimple(PARTICLETYPE);
    ParticleSimple(Particle*);
    ParticleSimple(ParticleSimple*);
    // ParticleSimple(const Particle*);
    // ParticleSimple(const ParticleSimple&);
    // ParticleSimple(const ParticleSimple&);

    // copy assignment
    // ParticleSimple& operator=(const ParticleSimple&&);

    bool operator==(Particle*) const;
    bool operator==(const Particle*) const;
    bool operator==(const ParticleSimple&) const;

    // members
    cartesian position;
    cartesian orientation;
    cartesian velocity;
    cartesian force;
    float mass;
    PARTICLETYPE type;
    unsigned int ID;
    enhance::ConcurrentDeque<enhance::observer_ptr<ParticleSimple>> regionQuery;
    enhance::observer_ptr<Particle> parent;
    tbb::atomic<bool> visited {false};
};