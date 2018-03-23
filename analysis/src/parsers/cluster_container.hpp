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

#include "cluster.hpp"



class ClusterParser
    : public Box<PERIODIC::ON>
    , public virtual ParameterDependentComponent
    , public enhance::ConcurrentDeque<Cluster>
{
public:
    typedef Cluster Cluster_t;
    typedef Cluster_t::type type;
    typedef Cluster_t::RegionQuery_t RegionQuery_t;

    template<typename iterator_t>
    void setTarget(iterator_t, iterator_t);

    template<typename FUNCTOR>
    void setupRegionQueries(FUNCTOR&&);

    void scan();

    float getOrder() const;
    Cluster& getLargest();
    const Cluster& getLargest() const;
    
    template<typename FUNCTOR> 
    enhance::observer_ptr<Cluster> getLargest(FUNCTOR&&);

protected:
    void DBSCAN_recursive(std::size_t = 1, float = 1.4);
    void expandCluster(Cluster&, type&, std::size_t);

    // enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
    tbb::concurrent_vector<type> particles {};
};



template<typename iterator_t>
void ClusterParser::setTarget(iterator_t _begin, iterator_t _end)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    clear();
    particles.clear();
    tbb::parallel_for_each(_begin, _end, [&](const auto& particle)
    {
        particles.emplace_back(particle.get());
    });
    assert(particles.size() == (std::size_t)std::abs(std::distance(_begin,_end)));

    // vesLOG("DBSCAN analysing " << particles.size() << " particles")
}



template<typename FUNCTOR>
void ClusterParser::setupRegionQueries(FUNCTOR&& condition)
{
    tbb::parallel_for_each(std::begin(particles), std::end(particles),[&](type& particle1)
    {
        particle1.regionQuery.clear();
        particle1.visited.store(false);
        std::for_each(std::begin(particles), std::end(particles),[&](type& particle2)
        {
            if(condition(particle1,particle2))
            {
                particle1.regionQuery.add(std::ref(particle2));
            }
        });
    });
}



template<typename FUNCTOR> 
enhance::observer_ptr<Cluster> ClusterParser::getLargest(FUNCTOR&& condition)
{
    enhance::observer_ptr<Cluster> largest_with_condition {nullptr};
    for(auto it = begin(); it != end(); ++it)
    {
        Cluster& cluster = *it;
        if(!largest_with_condition && condition(cluster))
            largest_with_condition = enhance::make_observer<Cluster>(std::addressof(cluster));
        else if(largest_with_condition && condition(cluster))
        {
            if(cluster.size() > largest_with_condition->size())
                largest_with_condition = enhance::make_observer<Cluster>(std::addressof(cluster));
        }
    }
    return largest_with_condition;
}