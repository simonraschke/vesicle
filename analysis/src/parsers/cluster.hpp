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
#include "definitions.hpp"
#include "surface_reconstruction.hpp"
#include <deque>
#include <ostream>
#include <tbb/parallel_for_each.h>



template<PERIODIC P = PERIODIC::ON>
class ClusterBase
    : public Box<P>
    , public virtual ParameterDependentComponent
    , public enhance::ConcurrentDeque<ParticleSimple>
{
public:
    typedef ParticleSimple type;
    typedef type::cartesian cartesian;
    typedef decltype(type::regionQuery)::member_t RegionQuery_t;

    // comparison operators
    template<PERIODIC P2> inline bool operator< (const ClusterBase<P2>& comp) const noexcept { return size() <  comp.size(); }
    template<PERIODIC P2> inline bool operator> (const ClusterBase<P2>& comp) const noexcept { return size() >  comp.size(); }
    
    template<typename FUNCTOR>
    void setupRegionQueries(FUNCTOR&&);

    template<PARTICLETYPE T>
    bool containsParticleType() const;

    template<PARTICLETYPE T>
    std::size_t numParticlesOfType() const;

    ClusterStructureParser& getStructure();
    const ClusterStructureParser& getStructure() const;

    cartesian getCenter() const;

    template<PERIODIC P2>
    friend std::ostream& operator<<(std::ostream&, const ClusterBase<P2>&);

    virtual ~ClusterBase() = default;

protected:
    ClusterBase() = default;

    std::unique_ptr<ClusterStructureParser> structure {nullptr};
};



class Subcluster
    : public ClusterBase<PERIODIC::OFF>
{

};



class Cluster
    : public ClusterBase<PERIODIC::OFF>
{
public:
    void scan();
    void rearrangeSubclusters();

    float getOrder() const;

    enhance::ConcurrentDeque<Subcluster> subclusters {};

protected:
    void DBSCAN_recursive(std::size_t, float);
    void expandCluster(Subcluster&, type&, std::size_t);
    void rearrangeSystemCopy();
    void rearrangeEdjucated();
    void rearrangeTrialAndError();
};



template<PERIODIC P>
template<typename FUNCTOR>
void ClusterBase<P>::setupRegionQueries(FUNCTOR&& condition)
{
    assert(size() % 8 == 0);

    tbb::parallel_for_each(begin(), end(),[&](type& particle1)
    {
        particle1.regionQuery.clear();
        particle1.visited.store(false);
        vesDEBUG("particle "<< particle1.ID << ": "<< particle1.position.format(ROWFORMAT))
        std::for_each(begin(), end(),[&](type& particle2)
        {
            if(condition(particle1,particle2))
            {
                vesDEBUG("compared successfully to "<< particle2.ID << ": " << particle2.position.format(ROWFORMAT));
                particle1.regionQuery.add(std::ref(particle2));
            }
        });
        if(particle1.regionQuery.size() > 12*getParameters().analysis_cluster_distance_threshold)
            throw std::logic_error("particle region query too large: "+std::to_string(particle1.regionQuery.size()));
    });
}



template<PERIODIC P>
template<PARTICLETYPE T>
bool ClusterBase<P>::containsParticleType() const
{
    return std::any_of(begin(), end(), [](const auto& p){ return T == p.type; });
}



template<PERIODIC P>
template<PARTICLETYPE T>
std::size_t ClusterBase<P>::numParticlesOfType() const
{
    return std::accumulate(begin(), end(), std::size_t(0), [](std::size_t i, const auto& p){ return p.type == T ? i+1 : i; });
}



template<PERIODIC P>
ClusterStructureParser& ClusterBase<P>::getStructure()
{
    if(!structure)
    {
        structure = std::make_unique<ClusterStructureParser>();
        structure->setParameters(getParameters());
        structure->setTarget(*this);
        structure->parse();
    }
    return *structure;
}



template<PERIODIC P>
const ClusterStructureParser& ClusterBase<P>::getStructure() const
{
    if(!structure)
    {
        throw std::logic_error(std::string("const calling ") + __func__ + ", but no structure was parsed yet");
    }
    return *structure;
}



template<PERIODIC P>
typename ClusterBase<P>::cartesian ClusterBase<P>::getCenter() const
{
    return std::accumulate(begin(), end(), cartesian(cartesian::Zero()), [&](const cartesian& c, const ParticleSimple& p){ return cartesian(c + p.position); }) / size();
}



template<PERIODIC P2>
std::ostream& operator<<(std::ostream& os, const ClusterBase<P2>& cluster)
{   
    os << "cluster of size " << cluster.size() << ": ";
    for(const auto& member : cluster)
        os << " " << member.ID;
    return os;
}