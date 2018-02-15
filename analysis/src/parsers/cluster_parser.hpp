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



template<PERIODIC P = PERIODIC::ON>
class ClusterParser
    : public Box<P>
    , public virtual ParameterDependentComponent
{
public:
    typedef ParticleSimple Particle_t;
    typedef enhance::observer_ptr<Particle_t> Particle_ptr_t;
    typedef enhance::ConcurrentDeque<Particle_ptr_t> Cluster_t;
    typedef std::deque<Cluster_t> ClusterList;

    // void setTarget(PARTICLERANGE*);
    template<typename iterator_t>
    void setTarget(iterator_t, iterator_t);

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



template<PERIODIC P>
template<typename iterator_t>
void ClusterParser<P>::setTarget(iterator_t _begin, iterator_t _end)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    particles.clear();
    clusters.clear();
    tbb::parallel_for_each(_begin, _end, [&](const auto& particle)
    {
        particles.emplace_back(particle.get());
    });

    vesLOG("DBSCAN analysing " << particles.size() << " particles")
}



template<PERIODIC P>
template<typename FUNCTOR>
void ClusterParser<P>::setupRegionQueries(FUNCTOR&& condition)
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


template<PERIODIC P>
void ClusterParser<P>::DBSCANrecursive(std::size_t min_points, float distance_threshold)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(range);

    // construct the regionQueries
    // this is done beforehand for performance reasons
    {
        const float sq_distance_threshold = distance_threshold*distance_threshold;
        setupRegionQueries([&](Particle_t& p1, Particle_t& p2)
        {
            return p1 == p2 || this->template squared_distance(p1.position, p2.position) <= sq_distance_threshold;
        });

    }

    // now DBSCAN
    {
        std::for_each(std::begin(particles), std::end(particles),[&](Particle_t& particle)
        {
            // only apporach particles, which were not scanned beforehand
            if(particle.visited.load() == false)
            {
                particle.visited.store(true);
            }
            else return;

            // add cluster if regionQuery is big enough
            if( particle.regionQuery.size() >= min_points )
            {
                clusters.emplace_back();
                auto ptr = enhance::make_observer<Particle_t>(&particle);
                clusters.back().add_if_new(ptr);

                // expanc cluster with regionQueries
                std::for_each(std::begin(particle.regionQuery), std::end(particle.regionQuery),[&](const Particle_ptr_t& other)
                {
                    if(ptr == other)
                        return;
                    
                    expandCluster(clusters.back(), other, min_points);
                });
            }
        });
    }

    vesLOG("DBSCAN found " << clusters.size() << " clusters with " << numParticles() << " particles")
}



template<PERIODIC P>
void ClusterParser<P>::expandCluster(Cluster_t& cluster, const Particle_ptr_t& particle, std::size_t min_points)
{
    // if(particle->visited.load() == false)
    // {
        particle->visited.store(true);
        cluster.add_if_new(particle);
    // }
    // else return;

    std::for_each(std::begin(particle->regionQuery), std::end(particle->regionQuery),[&](const Particle_ptr_t& other)
    {
        if(particle == other)
            return;

        if( other->visited.load() == false )
        {
            other->visited.store(true);

            // expand cluster if regionQuery is big enough
            if( particle->regionQuery.size() >= min_points )
            {
                // expand cluster with regionQueries
                std::for_each(std::begin(particle->regionQuery), std::end(particle->regionQuery),[&](const Particle_ptr_t& region_particle)
                {
                    expandCluster(cluster, region_particle, min_points);
                });
            }
        }
    });
}



template<PERIODIC P>
std::size_t ClusterParser<P>::numParticles() const
{
    return PARALLEL_REDUCE(std::size_t, clusters, [&](std::size_t i, const Cluster_t& cluster)
    {
        return i + cluster.size();
    });
}



template<PERIODIC P>
std::size_t ClusterParser<P>::numClusters() const
{
    return clusters.size();
}



template<PERIODIC P>
ClusterParser<P>::ClusterList::iterator ClusterParser<P>::begin()
{
    return std::begin(clusters);
}



template<PERIODIC P>
ClusterParser<P>::ClusterList::iterator ClusterParser<P>::end()
{
    return std::end(clusters);
}



template<PERIODIC P>
ClusterParser<P>::ClusterList::const_iterator ClusterParser<P>::begin() const
{
    return std::begin(clusters);
}



template<PERIODIC P>
ClusterParser<P>::ClusterList::const_iterator ClusterParser<P>::end() const
{
    return std::end(clusters);
}