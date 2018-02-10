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
#include "systems/box.hpp"
#include "vesicleIO/parameters.hpp"
#include <deque>
#include <tbb/parallel_for_each.h>


class ClusterContainer
    : public Box<PERIODIC::ON>
    , public virtual ParameterDependentComponent
{
public:
    void setTarget(PARTICLERANGE*);
    void setup(float = 1.4);
    
protected:

private:
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
    std::deque<Cluster> clusters;
    tbb::concurrent_vector<ParticleSimple> particles {};
};
