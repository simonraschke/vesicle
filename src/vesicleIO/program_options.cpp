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

#include "program_options.hpp"



// void ProgramOptions::read(int argc, const char* argv[])
// {
//     namespace po = boost::program_options;

//     optionsMap.clear();
//     std::string config_file_name;
    
//     po::options_description generalOptions("General Options");
//     generalOptions.add_options()
//         ("config", po::value<std::string>(&config_file_name), "read from config file")
//         ("help,h", "show help")
//         ("general.algorithm",  po::value<std::string>(), "[verlet,shakeVerlet,langevin,montecarlo]")
//         ("general.acceptance",  po::value<std::string>(), "[metropolis]")
//         ("general.interaction",  po::value<std::string>(), "[lj,alj]")
//         ("general.thermostat",  po::value<std::string>(), "[andersen]")
//     ;
    
//     po::options_description systemOptions("System Options");
//     systemOptions.add_options()
//         ("system.mobile,m", po::value<std::size_t>(), "number of mobile particles")
//         ("system.density,c", po::value<float>(), "particle density")
//         ("system.box.x", po::value<float>(), "box edge x")
//         ("system.box.y", po::value<float>(), "box edge y")
//         ("system.box.z", po::value<float>(), "box edge z")
//         ("system.temperature,t", po::value<float>(), "temperature")
//         ("system.timestep", po::value<float>(), "time step")
//         ("system.kappa,k", po::value<float>(), "kappa")
//         ("system.gamma,g", po::value<float>(), "gamma angle")
//         ("system.stepwidth_coordinates", po::value<float>(), "coordinates stepwidth (MonteCarlo only)")
//         ("system.stepwidth_orientation", po::value<float>(), "orientation stepwidth (MonteCarlo only)")
//         ("system.cell_min_edge", po::value<float>(), "minimum edge length of cell")
//         ("system.max_cells_dim", po::value<std::size_t>(), "maximum cells per dimension")
//     ;
    
//     po::options_description outputOptions("Output Options");
//     outputOptions.add_options()
//         ("output.traj",  po::value<std::string>(), "[gro]")
//         ("output.skip",  po::value<std::size_t>(), "print every .. step")
//     ;
    
//     po::options_description inputOptions("Input Options");
//     inputOptions.add_options()
//         ("input.traj",  po::value<std::string>(), "[none,gro]")
//         ("input.path",  po::value<boost::filesystem::path>(), "full or rel path to trajectory file")
//         ("input.frames",  po::value<std::string>(), "regular expression for frames to read")
//     ;

//     po::options_description analysisOptions("Analysis Options");
//     analysisOptions.add_options()
//         ("analysis.path", po::value<boost::filesystem::path>()->default_value("data.h5"), "HDF5 file to store all information")
//         ("analysis.full", po::bool_switch(), "enable full analysis")
//         ("analysis.epot", po::bool_switch(), "enable potential energy analysis")
//         ("analysis.cluster", po::bool_switch(), "enable full cluster analysis")
//         ("analysis.cluster_algorithm", po::value<std::string>()->default_value("DBSCAN"), "[DBSCAN]")
//         ("analysis.cluster_volume", po::bool_switch(), "enable cluster volume analysis")
//         ("analysis.cluster_histogram", po::bool_switch(), "enable cluster size histograms")
//     ;
    
//     po::options_description allOptions;
//     allOptions.add(generalOptions).add(systemOptions).add(outputOptions).add(inputOptions).add(analysisOptions);

//     po::store(po::command_line_parser(argc,argv).options(allOptions).run(),optionsMap);
//     po::notify(optionsMap);

//     PATH config_file_full_path = boost::filesystem::current_path() / config_file_name;

//     if(optionsMap.count("help"))
//     {
//         std::clog << '\n' << generalOptions;
//         std::clog << '\n' << systemOptions;
//         std::clog << '\n' << outputOptions;
//         std::clog << '\n' << inputOptions;
//         std::clog << '\n' << analysisOptions;
//         std::exit(EXIT_SUCCESS);
//     }
//     else if(boost::filesystem::exists(config_file_full_path) && !config_file_name.empty())
//     {
//         try
//         {
//             read_from_file(allOptions,optionsMap);
//         }
//         catch(boost::bad_any_cast& e)
//         {
//             vesCRITICAL(e.what())
//             std::exit(EXIT_FAILURE);
//         }
//     }
// }



// void ProgramOptions::read_from_file(boost::program_options::options_description& desc, boost::program_options::variables_map& vm )
// {
//     namespace po = boost::program_options;

//     IFSTREAM INPUT( vm["config"].as<std::string>() );
    
//     vm.clear();
    
//     po::store(po::parse_config_file(INPUT,desc),vm);
//     po::notify(vm);
// }
