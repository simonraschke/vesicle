#include "cluster.hpp"



void Cluster::scan()
{
    if(getParameters().analysis_cluster_algorithm == "DBSCAN")
        DBSCAN_recursive(getParameters().analysis_cluster_minimum_size, getParameters().analysis_cluster_distance_threshold);
    else
        vesCRITICAL("cluster algorithm " << getParameters().analysis_cluster_algorithm << " unknown")
}



void Cluster::DBSCAN_recursive(std::size_t min_points, float distance_threshold)
{
    vesDEBUG(__PRETTY_FUNCTION__);
    subclusters.clear();
    // construct the regionQueries
    // this is done beforehand for performance reasons
    {
        const float sq_distance_threshold = distance_threshold*distance_threshold;
        setupRegionQueries([&](type& p1, type& p2)
        {
            return (p1.position - p2.position).squaredNorm() <= sq_distance_threshold;
        });
    }

    // now DBSCAN
    vesDEBUG("particles in cluster " << size());
    std::for_each(begin(), end(),[&](type& particle)
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
            vesDEBUG("particle.regionQuery.size() " << particle.regionQuery.size());
            subclusters.add();
            subclusters.back().setParameters(getParameters());
            subclusters.back().add_if_new(particle);

            // expanc cluster with regionQueries
            std::for_each(std::begin(particle.regionQuery), std::end(particle.regionQuery),[&](ParticleSimple& other)
            {                
                if(!other.visited.load())
                    expandCluster(subclusters.back(), other, min_points);
            });
        }
    });

    assert( size() == (std::size_t) std::count_if(begin(), end(), [](const auto& particle){ return particle.visited.load(); }) );
    // vesLOG( size() << " and " << std::accumulate(subclusters.begin(), subclusters.end(), std::size_t(0), [](auto i, const auto& sub){ return i + sub.size();} ) );
    assert( size() == std::accumulate(subclusters.begin(), subclusters.end(), std::size_t(0), [](auto i, const auto& sub){ return i + sub.size();} ) );
    // assert( size() == std::accumulate(subclusters.begin(), subclusters.end(), std::size_t(0), [](auto i, const auto& sub){ return i + sub.size();} ));
    // vesLOG("DBSCAN found " << clusters.size() << " clusters with " << numParticles() << " particles")

    // rearrangeSubclusters();
    // std::for_each(std::begin(subclusters), std::end(subclusters), [](Subcluster& subcluster)
    // {
    //     subcluster.getStructure();
    // });
}



void Cluster::expandCluster(Subcluster& subcluster, type& particle, std::size_t min_points)
{
    if(particle.visited.load() == false)
    {
        particle.visited.store(true);
        subcluster.add_if_new(particle);
    }
    else return;

    if(particle.regionQuery.size() >= min_points)
    {
        std::for_each(std::begin(particle.regionQuery), std::end(particle.regionQuery),[&](ParticleSimple& other)
        {
            if( other.visited.load() == false )
            {
                expandCluster(subcluster, other, min_points);
            }
        });
    }
}



void Cluster::rearrangeSubclusters()
{
    // if(std::max_element(subclusters.begin(), subclusters.end())->size() == size())
    //     return;
    if(subclusters.size() == 1)
    {
        vesDEBUG("no " << __func__ << " needed, subclusters.size() " << subclusters.size());
        return;
    }

    rearrangeSystemCopy();

    // while(subclusters.size() > 1)
    // {
    //     rearrangeEdjucated();
    // }
    getStructure();
    if(subclusters.size() != 1)
        throw std::runtime_error("subclusters.size() == "+std::to_string(subclusters.size()));
}



void Cluster::rearrangeSystemCopy()
{
    if(subclusters.size() == 1)
    {
        vesDEBUG("no " << __func__ << " needed, subclusters.size() " << subclusters.size());
        return;
    }

    const auto original_size = size();
    decltype(members) origin_copy { std::begin(members), std::end(members) };

    {
        decltype(members) x_direction_copy;
        std::copy( std::begin(members), std::end(members), std::back_inserter(x_direction_copy) );
        for(auto& copy : x_direction_copy)
        {
            copy.position(0) += getLengthX();
        }
        std::copy(x_direction_copy.begin(), x_direction_copy.end(), std::back_inserter(members));
    }
    assert(size() == original_size*2);
    assert(std::abs(members[0].position(0) - members[original_size].position(0)) - getLengthX() < 1e-4);

    {
        decltype(members) y_direction_copy;
        std::copy( std::begin(members), std::end(members), std::back_inserter(y_direction_copy) );
        for(auto& copy : y_direction_copy)
        {
            copy.position(1) += getLengthY();
        }
        std::copy(y_direction_copy.begin(), y_direction_copy.end(), std::back_inserter(members));
    }
    assert(size() == original_size*4);
    assert(std::abs(members[0].position(1) - members[original_size*2].position(1)) - getLengthY() < 1e-4);

    {
        decltype(members) z_direction_copy;
        std::copy( std::begin(members), std::end(members), std::back_inserter(z_direction_copy) );
        for(auto& copy : z_direction_copy)
        {
            copy.position(2) += getLengthZ();
        }
        std::copy(z_direction_copy.begin(), z_direction_copy.end(), std::back_inserter(members));
    }
    assert(size() == original_size*8);
    assert(std::abs(members[0].position(2) - members[original_size*4].position(2)) - getLengthZ() < 1e-4);

    scan();

    const auto largest = std::max_element(subclusters.begin(), subclusters.end());

    // vesLOG("DBSCAN after copy of original_size " << original_size << ", now size "<< size())
    // for(auto& subcluster : subclusters)
    //     vesLOG("subcluster: size " << subcluster.size()  << "   " << subcluster)

    if(original_size == largest->size())
    {
        members = { std::begin(*largest), std::end(*largest) };
        vesDEBUG("members == largest subcluster of size " << size());
        scan();
    }
    else
    {
        vesWARNING("original size was "<< original_size << "  largest cluster has size " << largest->size())
        members.clear();
        members = { std::begin(origin_copy), std::end(origin_copy) };
        vesDEBUG("members == original_copy of size " << size());
        rearrangeEdjucated();
    }
    assert(size() == original_size);
}



void Cluster::rearrangeEdjucated()
{
}



float Cluster::getOrder() const
{
    const auto center = getCenter();
    return std::accumulate(begin(), end(), float(0), [&center](float order, const ParticleSimple& p)
    { 
        return order + (p.position-center).normalized().dot(p.orientation.normalized());
    }) / size();
}