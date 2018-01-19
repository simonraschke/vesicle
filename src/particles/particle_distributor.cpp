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



bool Distributor::conflicting_placement(PARTICLERANGE* range, PARTICLERANGE::value_type& p1)
{
    std::atomic<bool> conflict {false};
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p2) 
    {
        assert(p1 && p2);
        if(p1==p2 || conflict.load()) return;
        if(squared_distance(*p1,*p2) <= 1.122f*1.122f) conflict.store(true);
    });
    return conflict.load();
}



void RandomDistributor::operator()(PARTICLERANGE* range)
{
    assert(range);
    tbb::parallel_for_each(range->begin(), range->end(), [&](auto& p) 
    {
        assert(p);
        p->setCoords(randomCoords());
        p->setOrientation(randomOrientation());
    });

    assert(range);
    for(auto& p : *range)
    {
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
}



RandomDistributor::cartesian RandomDistributor::randomCoords() const
{
    return cartesian
    (
        enhance::random<cartesian::Scalar>(0.f,getLengthX()),
        enhance::random<cartesian::Scalar>(0.f,getLengthY()),
        enhance::random<cartesian::Scalar>(0.f,getLengthZ())
    );
}



RandomDistributor::cartesian RandomDistributor::randomOrientation() const
{
    return Eigen::Vector3f::Random();
}




void GridDistributor::operator()(PARTICLERANGE* range)
{   
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



void TrajectoryDistributor::operator()(PARTICLERANGE* range)
{   
    assert(range);
    if(getParameters().in_traj == std::string("gro"))
    {
        // get the lines of the last frame in regex range
        TrajectoryReaderGro reader;
        reader.setParameters(getParameters());
        reader.setPath(getParameters().in_traj_path);
        reader.readAllFrames();
        auto frame = reader.getMatches(getParameters().in_frames).rbegin()->second;

        for(std::size_t i = 0; i < reader.numParticles(); /**/)
        {
            Particle& particle = *((*range)[i]);
            const std::string line = frame[i+2];
            auto tokens = TrajectoryReaderGro::particleLineTokens(line);
            vesDEBUG("up    " << line);
            
            if(reader.isAnisotropic())
            {
                const std::string line2 = frame[i+3];
                auto tokens2 = TrajectoryReaderGro::particleLineTokens(line);
                vesDEBUG("bottom" << line2);
                {
                    cartesian position;
                    position(0) = (boost::lexical_cast<float>(tokens["pos x"]) + boost::lexical_cast<float>(tokens["pos x"])) /2;
                    position(1) = (boost::lexical_cast<float>(tokens["pos y"]) + boost::lexical_cast<float>(tokens["pos y"])) /2;
                    position(2) = (boost::lexical_cast<float>(tokens["pos z"]) + boost::lexical_cast<float>(tokens["pos z"])) /2;
                    particle.setCoords(position);
                }
                {
                    cartesian orientation;
                    orientation(0) = ( boost::lexical_cast<float>(tokens["pos x"]) - particle.coords()(0) );
                    orientation(1) = ( boost::lexical_cast<float>(tokens["pos y"]) - particle.coords()(1) );
                    orientation(2) = ( boost::lexical_cast<float>(tokens["pos z"]) - particle.coords()(2) );
                    particle.setOrientation(orientation);
                }
                {
                    cartesian velocity;
                    velocity(0) = (boost::lexical_cast<float>(tokens["vel x"]) + boost::lexical_cast<float>(tokens["vel x"])) /2;
                    velocity(1) = (boost::lexical_cast<float>(tokens["vel y"]) + boost::lexical_cast<float>(tokens["vel y"])) /2;
                    velocity(2) = (boost::lexical_cast<float>(tokens["vel z"]) + boost::lexical_cast<float>(tokens["vel z"])) /2;
                    particle.setVelocity(velocity);
                }
                i += 2;
            }
            else
            {
                {
                    cartesian position;
                    position(0) = boost::lexical_cast<float>(tokens["pos x"]);
                    position(1) = boost::lexical_cast<float>(tokens["pos y"]);
                    position(2) = boost::lexical_cast<float>(tokens["pos z"]);
                    particle.setCoords(position);
                }
                {
                    cartesian velocity;
                    velocity(0) = boost::lexical_cast<float>(tokens["vel x"]);
                    velocity(1) = boost::lexical_cast<float>(tokens["vel y"]);
                    velocity(2) = boost::lexical_cast<float>(tokens["vel z"]);
                    particle.setVelocity(velocity);
                }
                ++i;
            }
        }

    }
    else
    {
        throw std::logic_error("particles cannot be distributed via " + getParameters().in_traj);
    }
}