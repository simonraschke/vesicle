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

#include "particle_distributor.hpp"



bool Distributor::conflicting_placement(PARTICLERANGE* range, const PARTICLERANGE::value_type& p1) const
{
    std::atomic<bool> conflict {false};
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p2) 
    {
        assert(p1 && p2);
        if(p1==p2 || conflict.load()) return;
        if(squared_distance(*p1,*p2) <= std::pow(getParameters().LJsigma*1.122f,2)) conflict.store(true);
    });
    return conflict.load();
}



Distributor::cartesian Distributor::randomCoords() const
{
    return cartesian
    (
        enhance::random<cartesian::Scalar>(0.f,getLengthX()),
        enhance::random<cartesian::Scalar>(0.f,getLengthY()),
        enhance::random<cartesian::Scalar>(0.f,getLengthZ())
    );
}



Distributor::cartesian Distributor::randomOrientation() const
{
    return Eigen::Vector3f::Random();
}



void RandomDistributor::operator()(PARTICLERANGE* range)
{
    vesLOG("distributing particles randomly")
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p) 
    {
        assert(p);
        p->setCoords(randomCoords());
        p->setOrientation(randomOrientation());
    });

    struct IncHelper {};
    assert(range);
    for(auto& p : *range)
    {
        vesDEBUG(enhance::IncrementalNumberGenerator<IncHelper>()())
        std::size_t try_counter = 0;
        while(conflicting_placement(range,p))
        {
            assert(p);
            p->setCoords(randomCoords());
            if( ++try_counter > 1e6 )
            {
                vesWARNING("placing particle exceeded 1M tries")
                throw std::runtime_error("particle placement not possible");
            }
        }
    }
#ifndef NDEBUG
    for(auto& p1 : *range)
    {
        for(auto& p2: *range)
        {   
            if(p1==p2) continue;
            const auto sq_norm = this->distanceVector(p1->coords(),p2->coords()).squaredNorm();
            if(sq_norm < 1.f)
                vesCRITICAL("squared norm very small " << sq_norm)
        }
    }
#endif
}




void GridDistributor::operator()(PARTICLERANGE* range)
{   
    vesLOG("distributing particles on grid")
    std::size_t maxX = std::floor(getLengthX()/1.122f);
    std::size_t maxY = std::floor(getLengthY()/1.122f);
    std::size_t maxZ = std::floor(getLengthZ()/1.122f);

    // make a grid and a uniform distribution of its size
    GridGeometry grid(maxX,maxY,maxZ);
    grid.scale(cartesian(1.122f,1.122f,1.122f));
    std::uniform_int_distribution<std::size_t> dist(0,range->size()-1);
    // std::cout << maxX << " " << maxY << " " << maxZ << " = " << grid.points.size() << std::endl;
    
    // first set all to random point
    assert(range);
    const auto random_point = grid.points[dist(enhance::RandomEngine.pseudo_engine)];
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p) 
    {
        assert(p);
        p->setCoords(random_point);
    });

    // std::cout << "grid from: " << grid.points[0].format(ROWFORMAT) << "  to: " << grid.points.back().format(ROWFORMAT) << std::endl;

    // if grid has the same size as range populate all points
    if(range->size() == grid.points.size())
    {
        for(std::size_t i = 0; i < range->size(); ++i)
        {
            assert((*range)[i]);
            (*range)[i]->setCoords(grid.points[i]);
        }
    }
    // if smaller distribute randomly
    else if(range->size() < grid.points.size())
    {
        for(auto& p : *range)
        {
            // std::cout << "try: " << &p << "coords " << p->coords().format(ROWFORMAT) << std::endl;
            while(conflicting_placement(range,p))
            {
                assert(p);
                p->setCoords(grid.points[dist(enhance::RandomEngine.pseudo_engine)]);
                // std::cout << "coords after try: " << p->coords().format(ROWFORMAT) << std::endl;
            }
        }
    }
    else
    {
        throw std::logic_error("GridDistributor: More particles than grid points");
    }
}



void OsmoticSystemDistributor::operator()(PARTICLERANGE* range)
{
    const std::size_t num_frame = std::count_if(std::begin(*range), std::end(*range), [](const auto& particle){ return particle->getType() == PARTICLETYPE::FRAME;});
    const std::size_t num_mobile = std::count_if(std::begin(*range), std::end(*range), [](const auto& particle){ return particle->getType() == PARTICLETYPE::MOBILE;});
    const std::size_t num_osmotic = std::count_if(std::begin(*range), std::end(*range), [](const auto& particle){ return particle->getType() == PARTICLETYPE::OSMOTIC;});
    vesLOG("distributing particles for an osmotic pressure analysis system with " << num_frame << " frames, " << num_mobile << " mobiles, " << num_osmotic << " osmotic particles");

    const float optimum_distance = enhance::nth_root<6>(getParameters().LJsigma*2);
    const float radius = optimum_distance/(2.0*std::sin(getParameters().gamma));
    const float volume = enhance::sphere_volume(radius);
    std::size_t osmotic_bulk = std::round( getParameters().density * ((getParameters().x * getParameters().y * getParameters().z) - volume) );
    const std::size_t osmotic_inside = std::round( getParameters().osmotic_density_inside * volume);
    // vesLOG(getParameters().density << " " <<  getParameters().x << " " << getParameters().osmotic_density_inside);
    if(osmotic_inside+osmotic_bulk != getParameters().osmotic)
    {
        // vesCRITICAL("osmotic_inside " << osmotic_inside << "  osmotic_bulk " << osmotic_bulk << "PLEASE CHANGE DENSITY SLIGHTLY!")
        vesWARNING("changed osmotic_bulk from " << osmotic_bulk << " to " << getParameters().osmotic - osmotic_inside);
        osmotic_bulk = getParameters().osmotic - osmotic_inside;
    }

    vesLOG("optimum_distance " << optimum_distance)
    vesLOG("radius " << radius)
    vesLOG(osmotic_inside << " osmotic particles in micelle and " << osmotic_bulk << " in bulk")
    
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p) 
    {
        assert(p);
        p->setCoords(randomCoords());
        p->setOrientation(randomOrientation());
    });

    SphereGeometry sphere(getCenter(), radius, num_frame+num_mobile+2);
    assert(sphere.points.size() == num_frame+num_mobile+2);

    std::size_t sphere_counter = 0;
    std::size_t osmotic_counter = 0;
    for(std::size_t i = 0; i < range->size(); ++i)
    {
        if( sphere_counter == 1 || sphere_counter == sphere.points.size()-2 )
            ++sphere_counter;

        Particle& particle = *range->at(i);
        switch(particle.getType())
        {
            case UNDEFINED : 
                throw std::logic_error("Got particle of undefined type"); 
                break;
            case FRAME : 
                particle.setCoords(sphere.points[sphere_counter]); 
                particle.setOrientation(sphere.points[sphere_counter]-sphere.origin); 
                ++sphere_counter;
                break;
            case MOBILE : 
                particle.setCoords(sphere.points[sphere_counter]); 
                particle.setOrientation(sphere.points[sphere_counter]-sphere.origin); 
                ++sphere_counter;
                break;
            case OSMOTIC :
                if(Eigen::Vector3f vector; osmotic_counter < osmotic_inside)
                {
                    std::size_t try_counter = 0;
                    while((particle.coords() - getCenter()).norm() > radius-getParameters().kappa*1.1 || conflicting_placement(range,range->at(i)))
                    {
                        particle.setCoords(randomCoords());
                        if( ++try_counter > 1e6 )
                        {
                            vesWARNING("placing particle exceeded 1M tries")
                            throw std::runtime_error("particle placement not possible");
                        }
                    }
                    vesLOG("placed particle after " << try_counter << " tries")
                    // do
                    // {
                    //     vector = decltype(vector)(enhance::random<float>(0,getParameters().x) , enhance::random<float>(0,getParameters().y) , enhance::random<float>(0,getParameters().z));
                    // } while((vector - getCenter()).norm() > radius-getParameters().kappa*1.1 && conflicting_placement(range,range->at(i)));
                    // particle.setCoords(vector);
                }
                else
                {
                    std::size_t try_counter = 0;
                    while((particle.coords() - getCenter()).norm() < radius-getParameters().kappa*1.1 || conflicting_placement(range,range->at(i)))
                    {
                        particle.setCoords(randomCoords());
                        // vesLOG(particle.coords().format(ROWFORMAT));
                        if( ++try_counter > 1e6 )
                        {
                            vesWARNING("placing particle exceeded 1M tries")
                            throw std::runtime_error("particle placement not possible");
                        }
                    }
                    // do
                    // {
                    //     vector = decltype(vector)(enhance::random<float>(0,getParameters().x) , enhance::random<float>(0,getParameters().y) , enhance::random<float>(0,getParameters().z));
                    // } while((vector - getCenter()).norm() < radius+getParameters().kappa*1.1 && conflicting_placement(range,range->at(i)));
                    // particle.setCoords(vector); 
                }
                ++osmotic_counter;
                particle.setOrientation(cartesian(0,0,0));
                break;
            default :
                throw std::logic_error("encountered default in switch statement");
        }
    }
    
    tbb::parallel_for_each(range->begin(), range->end(), [](const auto& p)
    {   
        const auto coords = p->coords();
        if(std::isnan(coords(0))) { throw std::logic_error("particle coord has nan"); }
        if(std::isnan(coords(1))) { throw std::logic_error("particle coord has nan"); }
        if(std::isnan(coords(2))) { throw std::logic_error("particle coord has nan"); }
    });
}



void TrajectoryDistributorGro::operator()(PARTICLERANGE* range)
{   
    vesLOG("distributing particles from " << getParameters().in_traj_path)
    assert(range);
    if(getParameters().in_traj == std::string("gro"))
    {
        // get the lines of the last frame in regex range
        TrajectoryReaderGro reader;
        reader.setParameters(getParameters());
        reader.setPath(getParameters().in_traj_path);
        reader.readAllFrames();
        auto frame = reader.getMatches(getParameters().in_frames).rbegin()->second;
        // unsigned short anisotropicFactor = reader.isAnisotropic() ? 2 : 1;
        // std::size_t num_particles_lines = reader.numParticles() * anisotropicFactor;
        const std::size_t num_particles = reader.numParticles();
        if(num_particles != range->size())
            vesWARNING("PARTICLERANGE != num_particles from" << getParameters().in_traj_path)

        // vesCRITICAL("reader.isAnisotropic() " << std::boolalpha << reader.isAnisotropic());
        for(std::size_t i = 0; i < num_particles ; ++i)
        {
            PARTICLERANGE::value_type::element_type& particle = *((*range)[i]);
            const std::string line = frame[i*2+2];
            const auto tokens = reader.particleLineTokens(line);
            vesDEBUG("up    " << line);
            
            if(reader.isAnisotropic())
            {
                const std::string line2 = frame[i*2+3];
                const auto tokens2 = reader.particleLineTokens(line2);
                vesDEBUG("bottom" << line2);
                setupAnisotropicParticle(tokens, tokens2, particle);
                // i += 2;
            }
            else
            {
                setupIsotropicParticle(tokens, particle);
                // ++i;
            }
        }
        reader.clearAllFrames();
    }
    else
    {
        throw std::logic_error("particles cannot be distributed via " + getParameters().in_traj);
    }
}



void TrajectoryDistributorGro::setupAnisotropicParticle(const tokens_type& tokens, const tokens_type& tokens2, PARTICLERANGE::value_type::element_type& particle)
{
    {
        cartesian position;
        cartesian orientation;
        position(0) = (boost::lexical_cast<float>(tokens.at("pos x")) + boost::lexical_cast<float>(tokens2.at("pos x"))) /2;
        position(1) = (boost::lexical_cast<float>(tokens.at("pos y")) + boost::lexical_cast<float>(tokens2.at("pos y"))) /2;
        position(2) = (boost::lexical_cast<float>(tokens.at("pos z")) + boost::lexical_cast<float>(tokens2.at("pos z"))) /2;
        orientation(0) = boost::lexical_cast<float>(tokens.at("pos x")) - position(0);
        orientation(1) = boost::lexical_cast<float>(tokens.at("pos y")) - position(1);
        orientation(2) = boost::lexical_cast<float>(tokens.at("pos z")) - position(2);
        particle.setCoords(scaleDown(position));
        particle.setOrientation(orientation);
        vesDEBUG("position: " << position.format(ROWFORMAT) <<"  scaled down: " << particle.coords().format(ROWFORMAT) )
        vesDEBUG("orientation: " << orientation.format(ROWFORMAT) << "  direclty from particle: " << particle.orientation().format(ROWFORMAT) )
    }
    {
        cartesian velocity;
        velocity(0) = (boost::lexical_cast<float>(tokens.at("vel x")) + boost::lexical_cast<float>(tokens2.at("vel x"))) /2;
        velocity(1) = (boost::lexical_cast<float>(tokens.at("vel y")) + boost::lexical_cast<float>(tokens2.at("vel y"))) /2;
        velocity(2) = (boost::lexical_cast<float>(tokens.at("vel z")) + boost::lexical_cast<float>(tokens2.at("vel z"))) /2;
        particle.setVelocity(velocity);
    }
}



void TrajectoryDistributorGro::setupIsotropicParticle(const tokens_type& tokens, PARTICLERANGE::value_type::element_type& particle)
{
    {
        cartesian position;
        position(0) = boost::lexical_cast<float>(tokens.at("pos x"));
        position(1) = boost::lexical_cast<float>(tokens.at("pos y"));
        position(2) = boost::lexical_cast<float>(tokens.at("pos z"));
        particle.setCoords(position);
    }
    {
        cartesian velocity;
        velocity(0) = boost::lexical_cast<float>(tokens.at("vel x"));
        velocity(1) = boost::lexical_cast<float>(tokens.at("vel y"));
        velocity(2) = boost::lexical_cast<float>(tokens.at("vel z"));
        particle.setVelocity(velocity);
    }
}



void FrameGuidedGridDistributor::operator()(PARTICLERANGE* range)
{
    // const float radius = std::sqrt(1.1027)*getParameters().LJsigma/(std::sin(getParameters().gamma)*2);
    const float radius = std::pow(getParameters().LJsigma,1.0/6.0)/(2.0*std::sin(getParameters().gamma));
    const float dist_x = getLengthX()/getParameters().frame_guides_grid_edge;
    const float dist_y = getLengthY()/getParameters().frame_guides_grid_edge;
    const float dist_z = getLengthZ()/getParameters().frame_guides_grid_edge;

    vesLOG("distributing " << std::pow(getParameters().frame_guides_grid_edge,3) << " frame guides (" << getParameters().guiding_elements_each << " guiding elements each) on sphere grid and " << getParameters().mobile << " mobile particles randomly")
    vesLOG("radius of spheres is " << radius)
    SphereGridGeometry spheregrid(getParameters().frame_guides_grid_edge, radius, getParameters().guiding_elements_each);
    spheregrid.scale(cartesian(dist_x, dist_y, dist_z));
    spheregrid.shift(cartesian(dist_x/2, dist_y/2, dist_z/2));

    auto it = range->begin();
    for(const auto& sphere : spheregrid.spheres)
    {
        for(const auto& point : sphere.points)
        {
            Particle& particle = *(it->get());
            if(particle.getType() != FRAME)
            {
                throw std::logic_error("Didn't get PARTICLETYPE::FRAME, where it should have been");
            }
            else
            {
                particle.setCoords(point); 
                particle.setOrientation(point-sphere.origin);
            }

            // next particle
            std::advance(it,1);
        }
    }

    RandomDistributor random_dist;
    random_dist.setParameters(getParameters());

    tbb::parallel_for_each(it, range->end(), [&](auto& p) 
    {
        assert(p);
        p->setCoords(random_dist.randomCoords());
        p->setOrientation(random_dist.randomOrientation());
    });

    for(auto mobile_it = it; mobile_it != range->end(); ++mobile_it)
    {

        std::unique_ptr<Particle>& particle = *mobile_it;
        assert(particle);
        std::size_t try_counter = 0;
        while(conflicting_placement(range,particle))
        {
            assert(particle);
            particle->setCoords(random_dist.randomCoords());
            if( ++try_counter > 1e6 )
            {
                vesWARNING("placing particle exceeded 1M tries")
                throw std::runtime_error("particle placement not possible");
            }
        }
    }
}