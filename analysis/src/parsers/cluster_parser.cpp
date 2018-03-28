#include "cluster_parser.hpp"


void ClusterParser::scan()
{
    if(getParameters().analysis_cluster_algorithm == "DBSCAN")
        DBSCAN_recursive(getParameters().analysis_cluster_minimum_size, getParameters().analysis_cluster_distance_threshold);
    else
        vesCRITICAL("cluster algorithm " << getParameters().analysis_cluster_algorithm << " unknown")
}



void ClusterParser::DBSCAN_recursive(std::size_t min_points, float distance_threshold)
{
    vesDEBUG(__PRETTY_FUNCTION__)

    // construct the regionQueries
    // this is done beforehand for performance reasons
    clear();
    {
        const float sq_distance_threshold = distance_threshold*distance_threshold;
        setupRegionQueries([&](type& p1, type& p2)
        {
            return squared_distance(p1.position, p2.position) <= sq_distance_threshold;
        });
    }

    // now DBSCAN
    std::for_each(std::begin(particles), std::end(particles),[&](type& particle)
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
            add();
            back().setParameters(getParameters());
            back().add_if_new(particle);

            // expanc cluster with regionQueries
            std::for_each(std::begin(particle.regionQuery), std::end(particle.regionQuery),[&](ParticleSimple& other)
            {
                if(!other.visited.load())
                    expandCluster(back(), other, min_points);
            });
        }
    });

    assert(particles.size() == (std::size_t) std::count_if(std::begin(particles), std::end(particles), [](const auto& particle){ return particle.visited.load(); }) );
    vesLOG("DBSCAN found " << std::accumulate(begin(), end(), std::size_t(0), [](auto i, const auto& cluster){ return i + cluster.size();} ) << " particles with from " << particles.size());
    assert(particles.size() == std::accumulate(begin(), end(), std::size_t(0), [](auto i, const auto& cluster){ return i + cluster.size();} ) );
}



void ClusterParser::expandCluster(Cluster& cluster, type& particle, std::size_t min_points)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(particle.visited.load() == false)
    {
        particle.visited.store(true);
        cluster.add_if_new(particle);
    }
    else return;

    if(particle.regionQuery.size() >= min_points)
    {
        std::for_each(std::begin(particle.regionQuery), std::end(particle.regionQuery),[&](ParticleSimple& other)
        {
            if( other.visited.load() == false )
            {
                expandCluster(cluster, other, min_points);
            }
        });
    }
}



Cluster& ClusterParser::getLargest()
{
    return *std::max_element(this->begin(), this->end());
}



const Cluster& ClusterParser::getLargest() const
{
    return *std::max_element(this->begin(), this->end());
}



float ClusterParser::getOrder() const
{
    tbb::atomic<std::size_t> num_particles {0};
    return std::accumulate(begin(), end(), float(0), [&](float order, const Cluster& cluster)
    { 
        num_particles += cluster.size();
        return order + cluster.getOrder()*cluster.size(); 
    }) / num_particles;
}