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

#include "algorithm.hpp"
#include "monte_carlo_utility/metropolis.hpp"
// #include "monte_carlo_utility/cell_container.hpp"
#include "monte_carlo_utility/cell_based_algorithm.hpp"
#include "monte_carlo_utility/stepwidth_alignment_unit.hpp"
#include "particles/particle.hpp"
#include "enhance/triangular_matrix.hpp"
#include <memory>



class MonteCarlo
    : public Algorithm
{
    CellContainer<Particle> cells {};

public: 
    virtual void setup() override;
    virtual void step(const unsigned long& = 1) override;
    virtual std::size_t getCurrentStep() const override;

    virtual cell_container_type& getCells() final;
    
    StepwidhtAlignmentUnit sw_position;
    StepwidhtAlignmentUnit sw_orientation;

    float getEnergyMatrixSum();
    
protected:
    void doMCmove(const cell_type&);
    void doCoordinatesMove(const cell_type&, Particle&);
    void doOrientationMove(const cell_type&, Particle&);

    // float potentialEnergy(const std::unique_ptr<Particle>&) const;
    // float potentialEnergyInRegion(const cell_type&, const Particle&) const;
    std::size_t current_step {0}
;
    std::unique_ptr<enhance::TriangularMatrix<float>> energy_work {};
    std::unique_ptr<enhance::TriangularMatrix<float>> energy_old {};
};