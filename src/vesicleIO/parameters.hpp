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

#include "enhance/math_utility.hpp"
#include "program_options.hpp"
#include "definitions.hpp"
#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    // GENERAL
    std::string algorithm {};
    std::string acceptance {};
    std::string interaction {};
    std::string thermostat {};

    // SYSTEM
    std::size_t mobile {};
    float dt {};
    float x {};
    float y {};
    float z {};
    float temperature {};
    float kappa {};
    float gamma {};
    float stepwidth_coordinates {};
    float stepwidth_orientation {};
    float cell_min_edge {};
    std::size_t max_cells_dim {};

    // OUTPUT
    std::string traj {};
    std::size_t traj_skip {};

    // functions
    void setup();
    ProgramOptions programOptions {};
};





struct ParameterDependentComponent
{
    void setParameters(Parameters);
    const Parameters& getParameters() const;
    virtual ~ParameterDependentComponent() = default;

protected:
    Parameters& mutableAccess();
    ParameterDependentComponent() = default;

private:
    std::unique_ptr<Parameters> parameters {nullptr};
};