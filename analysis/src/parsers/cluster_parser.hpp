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

#include "systems/box.hpp"
#include "vesicleIO/parameters.hpp"
#include "enhance/concurrent_container.hpp"
#include "particles/particle_simple.hpp"
#include <deque>
#include <tbb/parallel_for_each.h>


class ClusterParser
    : public Box<PERIODIC::ON>
    , public virtual ParameterDependentComponent
{
public:
    typedef ParticleSimple Particle_t;
    typedef enhance::observer_ptr<Particle_t> Particle_ptr_t;
    typedef enhance::ConcurrentDeque<Particle_ptr_t> Cluster_t;
    typedef std::deque<Cluster_t> ClusterList;

    void setTarget(PARTICLERANGE*);
    void DBSCANrecursive(std::size_t = 1, float = 1.4);

    std::size_t numParticles() const;
    std::size_t numClusters() const;

    // iteration
    typename ClusterList::iterator begin();
    typename ClusterList::iterator end();
    typename ClusterList::const_iterator begin() const;
    typename ClusterList::const_iterator end() const;
    
protected:
    template<typename FUNCTOR>
    void setupRegionQueries(FUNCTOR&&);

    void expandCluster(Cluster_t&, const Particle_ptr_t&, std::size_t);

private:
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
    ClusterList clusters;
    tbb::concurrent_vector<Particle_t> particles {};
};




template<typename FUNCTOR>
void ClusterParser::setupRegionQueries(FUNCTOR&& condition)
{
    assert(range);
    // distance_threshold *= distance_threshold;

    tbb::parallel_for_each(std::begin(particles), std::end(particles),[&](Particle_t& particle1)
    {
        std::for_each(std::begin(particles), std::end(particles),[&](Particle_t& particle2)
        {
            if(condition(particle1,particle2))
            {
                particle1.regionQuery.add_if_new(enhance::make_observer<Particle_t>(&particle2));
            }
        });
    });
}