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

#include "enhance/math_utility.hpp"
#include "enhance/output_utility.hpp"
#include "definitions.hpp"
#include <memory>
#include <exception>
#include <cassert>
#include <string>
#include <limits>
#include <regex>
#include <tbb/task_scheduler_init.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>




// the actual Parameters class
struct Parameters
{
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    // GENERAL
    std::size_t cpu_threads {};
    std::string algorithm {};
    std::string acceptance {};
    std::string interaction {};
    std::string thermostat {};

    // SYSTEM
    std::size_t mobile {};
    std::size_t guiding_elements_each {};
    std::size_t frame_guides_grid_edge {};
    std::size_t guiding_elements_plane {};
    std::size_t osmotic {};
    std::size_t num_all_particles {};
    float osmotic_density_inside {};
    float density {};
    float x {};
    float y {};
    float z {};
    float temperature {};
    float dt {};
    float time_max {};
    float kappa {};
    float gamma {};
    float LJepsilon {};
    float LJsigma {};
    float sw_position_min {};
    float sw_position_max {};
    float sw_position_target {};
    float sw_orientation_min {};
    float sw_orientation_max {};
    float sw_orientation_target {};
    float cell_min_edge {};
    std::size_t max_cells_dim {};

    // OUTPUT
    std::string out_traj {};
    std::size_t out_traj_skip {};

    // INPUT
    std::string in_traj {};
    boost::filesystem::path in_traj_path {};
    std::regex in_frames {};

    // local memberfunctions 
    void setup();
    void read(int, const char* []);

    // static member functions
    static void read_from_file(boost::program_options::options_description&, boost::program_options::variables_map&);
    
protected:
    // the actual options map
    boost::program_options::variables_map optionsMap {};
};




// A parameter containing base class to derive from
// CARFUL: need to call setParameters on derived to not throw
// think carfully about virtual inheritance in multiple inheritance
struct ParameterDependentComponent
{
    // make_unique<Parameters> from copy of argument
    void setParameters(Parameters);

    // access underlying Parameters class
    // will throw if parameters private member is nullptr
    // prevent by calling setParameters before
    const Parameters& getParameters() const;

    // get a string map with relevant system parameters
    // to be used as attributes for hdf5 datasets
    std::map<std::string,std::string> systemAttributes() const;

    // destroy if derived is destroyed
    virtual ~ParameterDependentComponent() = default;

protected:
    // derived class is able to change values of parameters
    Parameters& mutableAccess();

    // only construct via inheritance
    ParameterDependentComponent() = default;

private:
    std::unique_ptr<Parameters> parameters {nullptr};
};
