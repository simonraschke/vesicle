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

#include "parameters.hpp"



void Parameters::read(int argc, const char* argv[])
{
    namespace po = boost::program_options;

    optionsMap.clear();
    std::string config_file_name;
    
    po::options_description generalOptions("General Options");
    generalOptions.add_options()
        ("config", po::value<std::string>(&config_file_name), "read from config file")
        ("help,h", "show help")
        ("general.algorithm",  po::value<std::string>(&algorithm)->default_value("montecarlo"), "[verlet,shakeVerlet,langevin,montecarlo]")
        ("general.acceptance",  po::value<std::string>(&acceptance)->default_value("metropolis"), "[metropolis]")
        ("general.interaction",  po::value<std::string>(&interaction)->default_value("alj"), "[lj,alj]")
        ("general.thermostat",  po::value<std::string>(&thermostat)->default_value("none"), "[none,andersen]")
    ;
    
    po::options_description systemOptions("System Options");
    systemOptions.add_options()
        ("system.mobile,m", po::value<std::size_t>(), "number of mobile particles")
        ("system.density,c", po::value<float>(), "particle density")
        ("system.box.x", po::value<float>(), "box edge x")
        ("system.box.y", po::value<float>(), "box edge y")
        ("system.box.z", po::value<float>(), "box edge z")
        ("system.temperature,t", po::value<float>(), "temperature")
        ("system.timestep", po::value<float>(&dt)->default_value(0.01), "time step")
        ("system.kappa,k", po::value<float>(&kappa)->default_value(1.0), "kappa")
        ("system.gamma,g", po::value<float>(&gamma)->default_value(11.5), "gamma angle")
        ("system.stepwidth_coordinates", po::value<float>(&stepwidth_coordinates)->default_value(0.2), "coordinates stepwidth (MonteCarlo only)")
        ("system.stepwidth_orientation", po::value<float>(&stepwidth_orientation)->default_value(0.2), "orientation stepwidth (MonteCarlo only)")
        ("system.cell_min_edge", po::value<float>(&cell_min_edge), "minimum edge length of cell")
        ("system.max_cells_dim", po::value<std::size_t>(&max_cells_dim)->default_value(20), "maximum cells per dimension")
    ;
    
    po::options_description outputOptions("Output Options");
    outputOptions.add_options()
        ("output.traj",  po::value<std::string>(&out_traj)->default_value("gro"), "[gro]")
        ("output.skip",  po::value<std::size_t>(&out_traj_skip)->default_value(10000), "print every .. step")
    ;
    
    po::options_description inputOptions("Input Options");
    inputOptions.add_options()
        ("input.traj",  po::value<std::string>(&in_traj)->default_value("none"), "[none,gro]")
        ("input.path",  po::value<boost::filesystem::path>(&in_traj_path)->default_value("trajectory.gro"), "full or rel path to trajectory file")
        ("input.frames",  po::value<std::string>()->default_value("-1"), "regular expression for frames to read")
    ;

    po::options_description analysisOptions("Analysis Options");
    analysisOptions.add_options()
        ("analysis.input", po::value<boost::filesystem::path>(&analysis_input)->default_value("trajectory.gro"), "trajectory file to analyse")
        ("analysis.path", po::value<boost::filesystem::path>(&analysis_path)->default_value("data.h5"), "HDF5 file to store all information")
        ("analysis.overwrite", po::bool_switch(&analysis_overwrite)->default_value(false), "overwrite existing file")
        ("analysis.frames", po::value<std::string>()->default_value("^[0-9]*$"), "regular expression for frames to analyse")
        ("analysis.full", po::bool_switch(&analysis_full)->default_value(true), "enable full analysis")
        ("analysis.epot", po::bool_switch(&analysis_epot)->default_value(false), "enable potential energy analysis")
        ("analysis.cluster", po::bool_switch(&analysis_cluster)->default_value(false), "enable full cluster analysis")
        ("analysis.cluster_algorithm", po::value<std::string>(&analysis_cluster_algorithm)->default_value("DBSCAN"), "[DBSCAN]")
        ("analysis.cluster_volume", po::bool_switch(&analysis_cluster_volume)->default_value(false), "enable cluster volume analysis")
        ("analysis.cluster_histogram", po::bool_switch(&analysis_cluster_histogram)->default_value(false), "enable cluster size histograms")
    ;
    
    po::options_description allOptions;
    allOptions.add(generalOptions).add(systemOptions).add(outputOptions).add(inputOptions).add(analysisOptions);

    po::store(po::command_line_parser(argc,argv).options(allOptions).run(),optionsMap);
    po::notify(optionsMap);

    PATH config_file_full_path = boost::filesystem::current_path() / config_file_name;

    if(optionsMap.count("help"))
    {
        std::clog << '\n' << generalOptions;
        std::clog << '\n' << systemOptions;
        std::clog << '\n' << outputOptions;
        std::clog << '\n' << inputOptions;
        std::clog << '\n' << analysisOptions;
        std::exit(EXIT_SUCCESS);
    }
    else if(boost::filesystem::exists(config_file_full_path) && !config_file_name.empty())
    {
        try
        {
            read_from_file(allOptions,optionsMap);
        }
        catch(boost::bad_any_cast& e)
        {
            vesCRITICAL(e.what())
            std::exit(EXIT_FAILURE);
        }
    }
}



void Parameters::read_from_file(boost::program_options::options_description& desc, boost::program_options::variables_map& vm )
{
    namespace po = boost::program_options;

    IFSTREAM INPUT( vm["config"].as<std::string>() );
    
    vm.clear();
    
    po::store(po::parse_config_file(INPUT,desc),vm);
    po::notify(vm);
}



void Parameters::setup()
{
    if(
        algorithm != "verlet" &&
        algorithm != "shakeVerlet" &&
        algorithm != "langevin" &&
        algorithm != "montecarlo"
    )
    {
        vesWARNING("general.algorithm not defined, choosing montecarlo")
        algorithm = "montecarlo";
    }



    if(
        acceptance != "metropolis"
    )
    {
        vesWARNING("general.acceptance not defined, choosing metropolis")
        acceptance = "metropolis";
    }



    if(
        interaction != "lj" &&
        interaction != "alj"
    )
    {
        vesWARNING("general.interaction not defined, choosing alj")
        interaction = "alj";
    }



    if(
        thermostat != "andersen" &&
        thermostat != "none"
    )
    {
        vesWARNING("general.thermostat not defined, choosing none")
        thermostat = "none";
    }


    // system configuration
    {
        unsigned short counter = 0;
        if(optionsMap.count("system.mobile")) ++counter;
        if(optionsMap.count("system.density")) ++counter;
        if(optionsMap.count("system.box.x") 
        && optionsMap.count("system.box.y") 
        && optionsMap.count("system.box.z")) ++counter;

        if(counter != 2)
        {
            vesCRITICAL("define 2 of the 3: system.mobile, system.density, system.box.{_}!")
        }
        else if(optionsMap.count("system.mobile") && optionsMap.count("system.density"))
        {
            vesLOG("system.mobile and system.density were set. Assuming cubic box")
            mobile = optionsMap["system.mobile"].as<std::size_t>();
            x = std::cbrt(static_cast<float>(mobile)/optionsMap["system.density"].as<float>());
            y = std::cbrt(static_cast<float>(mobile)/optionsMap["system.density"].as<float>());
            z = std::cbrt(static_cast<float>(mobile)/optionsMap["system.density"].as<float>());
        }
        else if(optionsMap.count("system.box.x") && optionsMap.count("system.density"))
        {
            x = optionsMap["system.box.x"].as<float>();
            y = optionsMap["system.box.y"].as<float>();
            z = optionsMap["system.box.z"].as<float>();
            mobile = std::round( optionsMap["system.density"].as<float>() * x * y * z);
        }
        else if(optionsMap.count("system.box.x") && optionsMap.count("system.box.y") && optionsMap.count("system.box.z")  && optionsMap.count("system.mobile"))
        {
            x = optionsMap["system.box.x"].as<float>();
            y = optionsMap["system.box.y"].as<float>();
            z = optionsMap["system.box.z"].as<float>();
            mobile = optionsMap["system.mobile"].as<std::size_t>();
            density = mobile/(x*y*z);
        }
        else 
            vesCRITICAL("UNKNOWN ERROR: maybe box not fully defined")


        if(optionsMap.count("system.temperature"))
            temperature = optionsMap["system.temperature"].as<float>();
        else 
            vesCRITICAL("system.temperature not defined")

        
        if(algorithm == std::string("montecarlo"))
            dt = 1;

        // degree to radians
        gamma = enhance::deg_to_rad(gamma);

        if(optionsMap.count("system.cell_min_edge"))
            cell_min_edge = optionsMap["system.cell_min_edge"].as<float>();
        else 
        {
            //TODO implement LJ sigma lazy ass
            vesWARNING("system.cell_min_edge not defined, choosing 3 Lennard Jones sigma")   
            cell_min_edge = 3.f * 1.f;
        }
    }

    // output configuration
    {
        if(
            out_traj != "gro" &&
            out_traj != "none"
        )
        {
            vesWARNING("output.out_traj not defined, no trajectory output")   
            out_traj = "none";
        }
    }

    // input configuration
    {
        if(
            in_traj != "gro" &&
            in_traj != "none"
        )
        {
            vesWARNING("input.in_traj not defined, no trajectory input")   
            in_traj = "none";
        }

        // if trajectory type set to none (expecting to be valid trajectory)
        // and path not existing
        // and path not none
        if(in_traj != std::string("none") && !boost::filesystem::exists(in_traj_path) && in_traj_path != std::string("none"))
        {
            vesCRITICAL("input.path=" << in_traj_path << " input trajectory not found, setting to \"none\"")
            in_traj_path = "none";
        }


        in_frames = optionsMap["input.frames"].as<std::string>();
    }

    // analysis configuration
    {
        // if trajectory type set to none (expecting to be valid trajectory)
        // and path not existing
        // and path not none
        if(enhance::splitAtDelimiter(analysis_input.string(),".").back() != "gro" )
        {
            vesCRITICAL("analysis.input=" << analysis_path << "  not a valid trajectory file format")
        }
        
        if(enhance::splitAtDelimiter(analysis_path.string(),".").back() != "h5" )
        {
            vesCRITICAL("analysis.path=" << analysis_path << "  not a .h5 file format")
        }
        
        analysis_frames = optionsMap["analysis.frames"].as<std::string>();
        
        if(analysis_full)
        {
            analysis_epot = true;
            analysis_cluster = true;
            analysis_cluster_volume = true;
            analysis_cluster_histogram = true;
        }

        if(!analysis_cluster)
        {
            analysis_cluster_volume = false;
            analysis_cluster_histogram = false;
        }
    }

    {
        vesLOG(__PRETTY_FUNCTION__)
        
        vesLOG("VALUE OVERVIEW")
        vesLOG("general.algorithm            " << algorithm )
        vesLOG("general.acceptance           " << acceptance )
        vesLOG("general.interaction          " << interaction )
        vesLOG("general.thermostat           " << thermostat )
        vesLOG("system.mobile                " << mobile )
        vesLOG("system.density               " << density )
        vesLOG("system.box.x                 " << x )
        vesLOG("system.box.y                 " << y )
        vesLOG("system.box.z                 " << z )
        vesLOG("system.timestep              " << dt )
        vesLOG("system.temperature           " << temperature )
        vesLOG("system.kappa                 " << kappa )
        vesLOG("system.gamma                 " << gamma )
        vesLOG("system.stepwidth_coordinates " << stepwidth_coordinates )
        vesLOG("system.stepwidth_orientation " << stepwidth_orientation )
        vesLOG("system.cell_min_edge         " << cell_min_edge )
        vesLOG("system.max_cells_dim         " << max_cells_dim )
        vesLOG("output.out_traj              " << out_traj )
        vesLOG("output.skip                  " << out_traj_skip )
        vesLOG("input.traj                   " << in_traj )
        vesLOG("input.path                   " << in_traj_path )
        vesLOG("input.frames                 " << optionsMap["input.frames"].as<std::string>() )
        vesLOG("analysis.input               " << analysis_input )
        vesLOG("analysis.path                " << analysis_path )
        vesLOG("analysis.frames              " << optionsMap["analysis.frames"].as<std::string>() )
        vesLOG("analysis.full                " << std::boolalpha << analysis_full )
        vesLOG("analysis.cluster             " << std::boolalpha << analysis_cluster )
        vesLOG("analysis.cluster_algorithm   " << analysis_cluster_algorithm )
        vesLOG("analysis.cluster_volume      " << std::boolalpha << analysis_cluster_volume )
        vesLOG("analysis.cluster_histogram   " << std::boolalpha << analysis_cluster_histogram )
    }
}



void ParameterDependentComponent::setParameters(Parameters prms)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    parameters = std::make_unique<Parameters>(prms);
}



const Parameters& ParameterDependentComponent::getParameters() const
{
    if(!parameters)
    {
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}



Parameters& ParameterDependentComponent::mutableAccess()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!parameters)
    {
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}
