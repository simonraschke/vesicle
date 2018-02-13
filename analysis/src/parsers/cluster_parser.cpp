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

#include "cluster_parser.hpp"



void ClusterParser::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);
    particles.clear();
    clusters.clear();
    tbb::parallel_for_each(std::begin(*range), std::end(*range),[&](const PARTICLERANGE::value_type& particle)
    {
        particles.emplace_back(particle.get());
    });


    vesLOG("DBSCAN analysing " << particles.size() << " particles")
}



void ClusterParser::DBSCANrecursive(std::size_t min_points, float distance_threshold)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(range);

    // construct the regionQueries
    // this is done beforehand for performance reasons
    {
        const float sq_distance_threshold = distance_threshold*distance_threshold;
        setupRegionQueries([&](Particle_t& p1, Particle_t& p2)
        {
            return p1 == p2 || squared_distance(p1.position, p2.position) <= sq_distance_threshold;
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



void ClusterParser::expandCluster(Cluster_t& cluster, const Particle_ptr_t& particle, std::size_t min_points)
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



std::size_t ClusterParser::numParticles() const
{
    return PARALLEL_REDUCE(std::size_t, clusters, [&](std::size_t i, const Cluster_t& cluster)
    {
        return i + cluster.size();
    });
}



std::size_t ClusterParser::numClusters() const
{
    return clusters.size();
}



ClusterParser::ClusterList::iterator ClusterParser::begin()
{
    return std::begin(clusters);
}



ClusterParser::ClusterList::iterator ClusterParser::end()
{
    return std::end(clusters);
}



ClusterParser::ClusterList::const_iterator ClusterParser::begin() const
{
    return std::begin(clusters);
}



ClusterParser::ClusterList::const_iterator ClusterParser::end() const
{
    return std::end(clusters);
}