/*  
*   Copyright 2017 Simon Raschke
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
#include "monte_carlo_utility/acceptance_adapter.hpp"
#include <tbb/tbb.h>
#include <eigen3/Eigen/Sparse>



// Simulation algorithm base class to derive from
// is Parameterdependent component

class Algorithm
    : public ParameterDependentComponent
{
public:
    virtual void setup();

    // set Parameters
    void setTarget(PARTICLERANGE*);

    template<typename I>
    void setInteraction();
    const std::unique_ptr<Interaction>& getInteraction() const;

    // this should be nullptr if not MonteCarlo
    template<typename A>
    void setAcceptance();
    const std::unique_ptr<AcceptanceAdapter>& getAcceptance() const;

    // execute
    virtual void step(const unsigned long& = 1) = 0;

    // necessary
    virtual ~Algorithm() = default;

protected:
    Algorithm() = default;

    virtual void updateCoords();
    virtual void updateOrientations();
    virtual void updateForces();
    virtual void updateVelocities();

    std::unique_ptr<Interaction> interaction {nullptr};
    std::unique_ptr<AcceptanceAdapter> acceptance {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};

    // std::unique_ptr<Eigen::SparseMatrix<float>> energy_matrix_old  {std::make_unique<Eigen::SparseMatrix<float>>()};
    // std::unique_ptr<Eigen::SparseMatrix<float>> energy_matrix_work {std::make_unique<Eigen::SparseMatrix<float>>()};

private:  

};



template<typename I>
void Algorithm::setInteraction()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    interaction = std::make_unique<I>();
    interaction->setParameters(getParameters());
    interaction->setup();
}



template<typename A>
void Algorithm::setAcceptance()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    acceptance = std::make_unique<A>();
    acceptance->setParameters(getParameters());
    // acceptance->setup();
}