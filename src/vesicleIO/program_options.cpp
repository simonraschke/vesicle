#include "program_options.hpp"



void ProgramOptions::read(int argc, const char* argv[])
{
    namespace po = boost::program_options;

    optionsMap.clear();
    std::string config_file_name;
    
    po::options_description generalOptions("General Options");
    generalOptions.add_options()
        ("config", po::value<std::string>(&config_file_name), "read from config file")
        ("help,h", "show help")
        ("general.algorithm",  po::value<std::string>(), "[verlet,shakeVerlet,langevin,montecarlo]")
        ("general.acceptance",  po::value<std::string>(), "[metropolis]")
        ("general.interaction",  po::value<std::string>(), "[lj,alj]")
        ("general.thermostat",  po::value<std::string>(), "[andersen]");
    
    po::options_description systemOptions("System Options");
    systemOptions.add_options()
        ("system.mobile,m", po::value<std::size_t>(), "number of mobile particles")
        ("system.density,c", po::value<float>(), "box edges x,y,z")
        ("system.box,b", po::value<std::vector<float>>(), "box edges x,y,z")
        ("system.temperature,t", po::value<float>(), "temperature")
        ("system.timestep", po::value<float>(), "time step")
        ("system.kappa,k", po::value<float>(), "kappa")
        ("system.gamma,g", po::value<float>(), "gamma angle")
        ("system.stepwidth_coordinates", po::value<float>(), "coordinates stepwidth (MonteCarlo only)")
        ("system.stepwidth_orientation", po::value<float>(), "orientation stepwidth (MonteCarlo only)");
    
    po::options_description outputOptions("Output Options");
    outputOptions.add_options()
        ("output.traj",  po::value<std::string>(), "[gro]")
        ("output.skip",  po::value<std::size_t>(), "print every .. step");

    
    po::options_description allOptions;
    allOptions.add(generalOptions).add(systemOptions).add(outputOptions);

    po::store(po::command_line_parser(argc,argv).options(allOptions).run(),optionsMap);
    po::notify(optionsMap);

    PATH config_file_full_path = boost::filesystem::current_path() / config_file_name;

    if(optionsMap.count("help"))
    {
        std::clog << '\n' << generalOptions;
        std::clog << '\n' << systemOptions;
        std::clog << '\n' << outputOptions;
        std::exit(EXIT_SUCCESS);
    }
    else if(boost::filesystem::exists(config_file_full_path) && !config_file_name.empty())
    {
        try
        {
            read_from_file(allOptions,optionsMap);
        }
        catch(boost::bad_any_cast e)
        {
            vesCRITICAL(e.what())
            std::exit(EXIT_FAILURE);
        }
    }
}



void ProgramOptions::read_from_file(boost::program_options::options_description& desc, boost::program_options::variables_map& vm )
{
    namespace po = boost::program_options;

    IFSTREAM INPUT( vm["config"].as<std::string>() );
    
    vm.clear();
    
    po::store(po::parse_config_file(INPUT,desc),vm);
    po::notify(vm);
}