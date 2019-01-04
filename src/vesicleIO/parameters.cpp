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
        ("general.cpu_threads",  po::value<std::size_t>(&cpu_threads)->default_value(tbb::task_scheduler_init::default_num_threads()), "maximum number of threads")
        ("general.algorithm",  po::value<std::string>(&algorithm)->default_value("montecarlo"), "[verlet,shakeVerlet,langevin,montecarlo]")
        ("general.acceptance",  po::value<std::string>(&acceptance)->default_value("metropolis"), "[metropolis]")
        ("general.interaction",  po::value<std::string>(&interaction)->default_value("alj"), "[lj,alj]")
        ("general.thermostat",  po::value<std::string>(&thermostat)->default_value("none"), "[none,andersen]")
        ("general.grand_canonical", po::value<bool>(&grand_canonical)->default_value(false), "flag for adding a grand canonical step to the algorithm")
    ;
    
    po::options_description systemOptions("System Options");
    systemOptions.add_options()
        ("system.mobile,m", po::value<std::size_t>(), "number of mobile particles")
        ("system.guiding_elements_each", po::value<std::size_t>(), "number of guiding elements per frame guide")
        ("system.frame_guides_grid_edge", po::value<std::size_t>(), "number of frame guides per dimension")
        ("system.guiding_elements_plane", po::value<bool>(&guiding_elements_plane)->default_value(false), "flag in order to make planar frame guide")
        ("system.plane_edge", po::value<float>(), "number of guiding elements in plane")
        // ("system.osmotic", po::value<std::size_t>(), "number of osmotic particles")
        ("system.osmotic_density_inside", po::value<float>(&osmotic_density_inside), "density of osmotic particles in bulk")
        ("system.density,c", po::value<float>(), "particle density")
        ("system.box.x", po::value<float>(), "box edge x")
        ("system.box.y", po::value<float>(), "box edge y")
        ("system.box.z", po::value<float>(), "box edge z")
        ("system.temperature,t", po::value<float>(), "temperature")
        ("system.timestep", po::value<float>(&dt)->default_value(0.01), "time step")
        ("system.time_max", po::value<float>(&time_max)->default_value(FLT_MAX), "time step")
        ("system.ljepsilon,k", po::value<float>(&LJepsilon)->default_value(1.0), "kappa")
        ("system.ljsigma,k", po::value<float>(&LJsigma)->default_value(1.0), "kappa")
        ("system.kappa,k", po::value<float>(&kappa)->default_value(1.0), "kappa")
        ("system.gamma,g", po::value<float>(&gamma)->default_value(11.5), "gamma angle")
        ("system.sw_position_min", po::value<float>(&sw_position_min)->default_value(0.05), "position stepwidth min (MonteCarlo only)")
        ("system.sw_position_max", po::value<float>(&sw_position_max)->default_value(1.0), "position stepwidth min (MonteCarlo only)")
        ("system.sw_position_target", po::value<float>(&sw_position_target)->default_value(0.3), "position stepwidth min (MonteCarlo only)")
        ("system.sw_orientation_min", po::value<float>(&sw_orientation_min)->default_value(M_PI_4/20), "orientation stepwidth min (MonteCarlo only)")
        ("system.sw_orientation_max", po::value<float>(&sw_orientation_max)->default_value(M_PI_4), "orientation stepwidth max (MonteCarlo only)")
        ("system.sw_orientation_target", po::value<float>(&sw_orientation_target)->default_value(0.3), "orientation stepwidth max (MonteCarlo only)")
        ("system.cell_min_edge", po::value<float>(&cell_min_edge), "minimum edge length of cell")
        ("system.max_cells_dim", po::value<std::size_t>(&max_cells_dim)->default_value(20), "maximum cells per dimension")
    ;
    
    po::options_description outputOptions("Output Options");
    outputOptions.add_options()
        ("output.traj",  po::value<std::string>(&out_traj)->default_value("h5"), "[gro,h5]")
        ("output.skip",  po::value<std::size_t>(&out_traj_skip)->default_value(10000), "print every .. steps")
    ;
    
    po::options_description inputOptions("Input Options");
    inputOptions.add_options()
        ("input.traj",  po::value<std::string>(&in_traj)->default_value("none"), "[none,gro]")
        ("input.path",  po::value<boost::filesystem::path>(&in_traj_path)->default_value("trajectory.gro"), "full or rel path to trajectory file")
        ("input.frames",  po::value<std::string>()->default_value("-1"), "regular expression for frames to read")
    ;
    
    po::options_description allOptions;
    allOptions.add(generalOptions).add(systemOptions).add(outputOptions).add(inputOptions);

    po::store(po::command_line_parser(argc,argv).options(allOptions).run(),optionsMap);
    po::notify(optionsMap);

    PATH config_file_full_path = boost::filesystem::current_path() / config_file_name;

    if(optionsMap.count("help"))
    {
        std::clog << '\n' << generalOptions;
        std::clog << '\n' << systemOptions;
        std::clog << '\n' << outputOptions;
        std::clog << '\n' << inputOptions;
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
        acceptance == "metropolis"
    )
    {
        // yeah i know
    }
    else
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

        guiding_elements_each = optionsMap.count("system.guiding_elements_each") ? optionsMap["system.guiding_elements_each"].as<std::size_t>() : 0;
        frame_guides_grid_edge = optionsMap.count("system.frame_guides_grid_edge") ? optionsMap["system.frame_guides_grid_edge"].as<std::size_t>() : 0;
        guiding_elements_each = optionsMap.count("system.guiding_elements_plane") ? optionsMap["system.guiding_elements_each"].as<std::size_t>() : 0;
        osmotic = optionsMap.count("system.osmotic") ? optionsMap["system.osmotic"].as<std::size_t>() : 0;
        if( std::abs(osmotic_density_inside) < 1e-6 )
        {
            osmotic = 0;
        }

        if(counter != 2)
        {
            vesCRITICAL("define 2 of the 3: system.mobile, system.density, system.box.{_}!")
        }
        else if(optionsMap.count("system.mobile") && optionsMap.count("system.density"))
        {
            vesLOG("system.mobile and system.density were set. Assuming cubic box")
            mobile = optionsMap["system.mobile"].as<std::size_t>();
            density = optionsMap["system.density"].as<float>();
            // calc first time without osmotic particles
            num_all_particles = guiding_elements_plane ? guiding_elements_each + mobile : mobile + guiding_elements_each * std::pow(frame_guides_grid_edge,3);
            x = std::cbrt(static_cast<float>(num_all_particles)/density);
            y = std::cbrt(static_cast<float>(num_all_particles)/density);
            z = std::cbrt(static_cast<float>(num_all_particles)/density);
            if(std::abs(osmotic_density_inside) > 1e-6)
            {
                vesLOG("calculating osmotic particles")
                // FIXME: this is bullcrap
                const float optimum_distance = enhance::nth_root<6>(LJsigma*2);
                // const float radius = optimum_distance/(2.0*std::sin(enhance::deg_to_rad(gamma)));
                const float radius = optimum_distance*std::sqrt(1.1027*num_all_particles)/4;
                x = 2.0*radius + 6.0*LJsigma;
                y = 2.0*radius + 6.0*LJsigma;
                z = 2.0*radius + 6.0*LJsigma;
                const float volume = enhance::sphere_volume(radius);
                osmotic = std::round( density * ((x * y * z) - volume) ) + std::round( osmotic_density_inside * volume); 
                // calc first time with osmotic particles
                num_all_particles = guiding_elements_plane ? guiding_elements_each + mobile + osmotic : mobile + osmotic + guiding_elements_each * std::pow(frame_guides_grid_edge,3);
            }
        }
        else if(optionsMap.count("system.box.x") && optionsMap.count("system.density"))
        {
            x = optionsMap["system.box.x"].as<float>();
            y = optionsMap["system.box.y"].as<float>();
            z = optionsMap["system.box.z"].as<float>();
            density = optionsMap["system.density"].as<float>();
            std::size_t num_all_particles_minus_mobile = guiding_elements_plane ? osmotic + guiding_elements_each : osmotic + guiding_elements_each * std::pow(frame_guides_grid_edge,3) + guiding_elements_plane;
            mobile = std::round( density * x * y * z) - num_all_particles_minus_mobile;
            const float optimum_distance = enhance::nth_root<6>(LJsigma*2);
            const float radius = optimum_distance/(2.0*std::sin(enhance::deg_to_rad(gamma)));
            const float volume = enhance::sphere_volume(radius);
            osmotic = std::round( density * ((x * y * z) - volume) ) + std::round( osmotic_density_inside * volume); 
            num_all_particles = guiding_elements_plane ? guiding_elements_each + mobile + osmotic: mobile + osmotic + guiding_elements_each * std::pow(frame_guides_grid_edge,3);;
        }
        else if(optionsMap.count("system.box.x") && optionsMap.count("system.box.y") && optionsMap.count("system.box.z")  && optionsMap.count("system.mobile"))
        {
            x = optionsMap["system.box.x"].as<float>();
            y = optionsMap["system.box.y"].as<float>();
            z = optionsMap["system.box.z"].as<float>();
            mobile = optionsMap["system.mobile"].as<std::size_t>();
            const float optimum_distance = enhance::nth_root<6>(LJsigma*2);
            const float radius = optimum_distance/(2.0*std::sin(enhance::deg_to_rad(gamma)));
            const float volume = enhance::sphere_volume(radius);
            osmotic = std::round( density * ((x * y * z) - volume) ) + std::round( osmotic_density_inside * volume); 
            num_all_particles = guiding_elements_plane ? guiding_elements_each + mobile + osmotic: mobile + osmotic + guiding_elements_each * std::pow(frame_guides_grid_edge,3);
            density = num_all_particles/(x*y*z);
        }
        else 
        {
            vesCRITICAL("UNKNOWN ERROR: maybe box not fully defined")
        }

        if(optionsMap.count("system.plane_edge"))
            plane_edge = optionsMap["system.plane_edge"].as<float>();

        // again set osmotic to 0 if no density inside micelle is 0
        if( std::abs(osmotic_density_inside) < 1e-6 )
        {
            osmotic = 0;
        }

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
            out_traj != "h5" &&
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
            out_traj != "h5" &&
            in_traj != "none"
        )
        {
            vesWARNING("input.in_traj not defined, no trajectory input")
            in_traj = "none";
        }

        if( in_traj == "none" )
            GLOBAL::getInstance().mode.store(GLOBAL::NEWRUN);
        else
            GLOBAL::getInstance().mode.store(GLOBAL::RESTART);

        // if trajectory type set to none (expecting to be valid trajectory)
        // and path not existing
        // and path not none
        if(in_traj != std::string("none") && !boost::filesystem::exists(in_traj_path) && in_traj_path != std::string("none"))
        {
            vesWARNING("input.path=" << in_traj_path << " input trajectory not found, setting to \"none\"")
            in_traj_path = "none";
            GLOBAL::getInstance().mode.store(GLOBAL::NEWRUN);
        }


        in_frames = optionsMap["input.frames"].as<std::string>();
    }

    {
        vesLOG(__PRETTY_FUNCTION__)
        
        vesLOG("VALUE OVERVIEW")
        vesLOG("general.cpu_threads                 " << cpu_threads )
        vesLOG("general.algorithm                   " << algorithm )
        vesLOG("general.acceptance                  " << acceptance )
        vesLOG("general.interaction                 " << interaction )
        vesLOG("general.thermostat                  " << thermostat )
        vesLOG("general.grand_canonical             " << grand_canonical )
        vesLOG("system.mobile                       " << mobile )
        vesLOG("system.frame_guides_grid_edge       " << frame_guides_grid_edge )
        vesLOG("system.guiding_elements_each        " << guiding_elements_each )
        vesLOG("system.guiding_elements_plane       " << guiding_elements_plane )
        vesLOG("system.plane_edge                   " << plane_edge )
        vesLOG("system.osmotic                      " << osmotic )
        vesLOG("system.osmotic_density_inside       " << osmotic_density_inside )
        vesLOG("system.num_all_particles            " << num_all_particles )
        vesLOG("system.density                      " << density )
        vesLOG("system.box.x                        " << x )
        vesLOG("system.box.y                        " << y )
        vesLOG("system.box.z                        " << z )
        vesLOG("system.timestep                     " << dt )
        vesLOG("system.time_max                     " << time_max )
        vesLOG("system.temperature                  " << temperature )
        vesLOG("system.kappa                        " << kappa )
        vesLOG("system.gamma                        " << gamma )
        vesLOG("system.ljepsilon                    " << LJepsilon )
        vesLOG("system.ljsigma                      " << LJsigma )
        vesLOG("system.sw_position_min              " << sw_position_min )
        vesLOG("system.sw_position_max              " << sw_position_max )
        vesLOG("system.sw_position_target           " << sw_position_target )
        vesLOG("system.sw_orientation_min           " << sw_orientation_min )
        vesLOG("system.sw_orientation_max           " << sw_orientation_max )
        vesLOG("system.sw_orientation_target        " << sw_orientation_target )
        vesLOG("system.cell_min_edge                " << cell_min_edge )
        vesLOG("system.max_cells_dim                " << max_cells_dim )
        vesLOG("output.out_traj                     " << out_traj )
        vesLOG("output.skip                         " << out_traj_skip )
        vesLOG("input.traj                          " << in_traj )
        vesLOG("input.path                          " << in_traj_path )
        vesLOG("input.frames                        " << optionsMap["input.frames"].as<std::string>() )
    }
}



std::map<std::string,std::string> ParameterDependentComponent::systemAttributes() const
{
    std::map<std::string,std::string> attributes;
    attributes["system.mobile"] = boost::lexical_cast<std::string>(getParameters().mobile);
    attributes["system.guiding_elements_each"] = boost::lexical_cast<std::string>(getParameters().guiding_elements_each);
    attributes["system.frame_guides_grid_edge"] = boost::lexical_cast<std::string>(getParameters().frame_guides_grid_edge);
    attributes["system.osmotic"] = boost::lexical_cast<std::string>(getParameters().osmotic);
    attributes["system.density"] = boost::lexical_cast<std::string>(getParameters().density);
    attributes["system.box.x"] = boost::lexical_cast<std::string>(getParameters().x);
    attributes["system.box.y"] = boost::lexical_cast<std::string>(getParameters().y);
    attributes["system.box.z"] = boost::lexical_cast<std::string>(getParameters().z);
    attributes["system.temperature"] = boost::lexical_cast<std::string>(getParameters().temperature);
    attributes["system.kappa"] = boost::lexical_cast<std::string>(getParameters().kappa);
    attributes["system.gamma"] = boost::lexical_cast<std::string>(getParameters().gamma);

    if(getParameters().interaction == "lj" || getParameters().interaction == "alj")
    {
        attributes["system.ljepsilon"] = boost::lexical_cast<std::string>(getParameters().LJepsilon);
        attributes["system.ljsigma"] = boost::lexical_cast<std::string>(getParameters().LJsigma);
    }
    return attributes;
}



void ParameterDependentComponent::setParameters(Parameters prms)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    parameters = std::make_unique<Parameters>(prms);
}


#include "enhance/stacktrace.cxx"
const Parameters& ParameterDependentComponent::getParameters() const
{
    if(!parameters)
    {
        #if defined(NDEBUG)
        vesLOG(Backtrace());
        #endif // NDEBUG
        
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}



Parameters& ParameterDependentComponent::mutableAccess()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!parameters)
    {
        #if defined(NDEBUG)
        vesLOG(Backtrace());
        #endif // NDEBUG
        
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}
