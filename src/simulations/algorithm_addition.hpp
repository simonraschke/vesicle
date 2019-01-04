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

#include "definitions.hpp"
#include "enhance/observer_ptr.hpp"
#include "vesicleIO/parameters.hpp"
#include "particles/particle.hpp"
#include "interactions/interaction.hpp"
// #include "monte_carlo_utility/acceptance_adapter.hpp"
#include "simulations/monte_carlo_utility/metropolis.hpp"
#include "particles/particle_distributor.hpp"
// #include "interactions/angular_lennard_jones.hpp"
// #include "interactions/lennard_jones.hpp"
#include <tbb/tbb.h>
#include <boost/progress.hpp>




class System;

class AlgorithmAddition
    : public ParameterDependentComponent
{
public:
    // void setTarget(PARTICLERANGE*);
    void setParentSystem(System*);
    virtual void apply() = 0;
    virtual void setup();

    // necessary
    virtual ~AlgorithmAddition() = default;

protected:
    AlgorithmAddition() = default;

    enhance::observer_ptr<System> parent_system {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};

private:  

};




class GrandCanonicalAddition
    : public AlgorithmAddition
{
public:
    virtual void apply() override;
    virtual void setup() override;

protected:
    bool tryRemove();
    bool tryInsert();
    bool conflicting_placement() const;
    
    MetropolisAcceptanceAdapter metropolis;
    std::unique_ptr<RandomDistributor> distributor {};
};



#include "systems/system.hpp"