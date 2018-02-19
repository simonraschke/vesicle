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

#include "vesicleIO/parameters.hpp"
#include "vesicleIO/gro_reader.hpp"
#include "translator/snapshot_translator.hpp"
#include "potential_energy.hpp"
#include "cluster_parser.hpp"
#include "surface_reconstruction.hpp"
#include "systems/controller.hpp"
#include "definitions.hpp"


#undef H5_USE_BOOST
#define H5_USE_BOOST
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include "H5File.hpp"


class DataCollector
    : public ParameterDependentComponent
{    
public:
    DataCollector();
    void setup();
    void collect();
    void write();

private:
    void try_potential_energy();
    void try_cluster();

    boost::filesystem::path working_dir;
    std::unique_ptr<HighFive::File> FILE {nullptr};
    std::unique_ptr<TrajectoryReader> reader {nullptr};
    AnisotropicSnapshotTranslatorGro translator;

    // collecting time_elapsed
    std::deque<float> timepoints {};

    // collecting potential energy
    std::deque<float> potential_energies {};

    // collecting cluster structures
    ClusterParser<PERIODIC::ON> clusters;
    std::deque<float> cluster_volumes {};
    std::deque<float> cluster_surface_areas {};
};