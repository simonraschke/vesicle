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

#include "cluster_container.hpp"



void ClusterContainer::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);
    particles.clear();
    tbb::parallel_for_each(std::begin(*range), std::end(*range),[&](const PARTICLERANGE::value_type& particle)
    {
        particles.emplace_back(particle.get());
    });
}



void ClusterContainer::setup(float distance_threshold)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(range);

    distance_threshold *= distance_threshold;
    tbb::parallel_for_each(std::begin(particles), std::end(particles),[&](ParticleSimple& particle1)
    {
        std::for_each(std::begin(particles), std::end(particles),[&](ParticleSimple& particle2)
        {
            if(particle1==particle2 || squared_distance(particle1.position, particle2.position)<=distance_threshold)
            {
                particle1.neighbourlist.emplace_back(particle2);
            }
        });
    });
}