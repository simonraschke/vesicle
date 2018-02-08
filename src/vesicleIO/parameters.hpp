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
#include "program_options.hpp"
#include "definitions.hpp"
#include <memory>
#include <exception>
#include <cassert>
#include <string>




// the actual Parameters class
struct Parameters
{
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    // GENERAL
    std::string algorithm {};
    std::string acceptance {};
    std::string interaction {};
    std::string thermostat {};

    // SYSTEM
    std::size_t mobile {};
    float density {};
    float x {};
    float y {};
    float z {};
    float temperature {};
    float dt {};
    float kappa {};
    float gamma {};
    float stepwidth_coordinates {};
    float stepwidth_orientation {};
    float cell_min_edge {};
    std::size_t max_cells_dim {};

    // OUTPUT
    std::string out_traj {};
    std::size_t out_traj_skip {};

    // INPUT
    std::string in_traj {};
    boost::filesystem::path in_traj_path {};
    std::regex in_frames {};

    // ANALYSIS
    boost::filesystem::path analysis_path {};
    bool analysis_overwrite {};
    std::regex analysis_frames {};
    bool analysis_full {};
    bool analysis_epot {};
    bool analysis_cluster {};
    std::string analysis_cluster_algorithm {};
    bool analysis_cluster_volume {};
    bool analysis_cluster_histogram {};

    // local memberfunctions 
    void setup();
    void read(int, const char* []);

    // static member functions
    static void read_from_file(boost::program_options::options_description&, boost::program_options::variables_map&);

    // TODO get rid of
    ProgramOptions programOptions {};
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