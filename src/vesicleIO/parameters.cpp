#include "parameters.hpp"



void Parameters::setup()
{
    if(programOptions.optionsMap.count("general.algorithm"))
        algorithm = programOptions.optionsMap["general.algorithm"].as<std::string>();
    else 
    {
        vesWARNING("general.algorithm not defined, choosing verlet")
        algorithm = "verlet";
    }

    if(programOptions.optionsMap.count("general.interaction"))
        interaction = programOptions.optionsMap["general.interaction"].as<std::string>();
    else 
    {
        vesWARNING("general.interaction not defined, choosing lj")
        interaction = "lj";
    }

    if(programOptions.optionsMap.count("general.thermostat"))
        thermostat = programOptions.optionsMap["general.thermostat"].as<std::string>();
    else 
    {
        vesWARNING("general.thermostat not defined, choosing andersen")
        thermostat = "andersen";
    }


    // system configuration
    {
        unsigned short counter = 0;
        if(programOptions.optionsMap.count("system.mobile")) ++counter;
        if(programOptions.optionsMap.count("system.density")) ++counter;
        if(programOptions.optionsMap.count("system.box")) ++counter;

        if(counter != 2)
            vesCRITICAL("define 2 of the 3: system.mobile, system.density, system.box!")
        else if(programOptions.optionsMap.count("system.mobile") && programOptions.optionsMap.count("system.density"))
        {
            vesLOG("system.mobile and system.density were set. Assuming cubic box")
            mobile = programOptions.optionsMap["system.mobile"].as<std::size_t>();
            x = std::cbrt(static_cast<float>(mobile)/programOptions.optionsMap["system.density"].as<float>());
            y = std::cbrt(static_cast<float>(mobile)/programOptions.optionsMap["system.density"].as<float>());
            z = std::cbrt(static_cast<float>(mobile)/programOptions.optionsMap["system.density"].as<float>());
        }
        else if(programOptions.optionsMap.count("system.box") && programOptions.optionsMap.count("system.density"))
        {
            x = programOptions.optionsMap["system.box"].as<std::vector<float>>()[0];
            y = programOptions.optionsMap["system.box"].as<std::vector<float>>()[1];
            z = programOptions.optionsMap["system.box"].as<std::vector<float>>()[2];
            mobile = std::round( programOptions.optionsMap["system.density"].as<float>() * x * y * z);
        }
        else if(programOptions.optionsMap.count("system.box") && programOptions.optionsMap.count("system.mobile"))
        {
            x = programOptions.optionsMap["system.box"].as<std::vector<float>>()[0];
            y = programOptions.optionsMap["system.box"].as<std::vector<float>>()[1];
            z = programOptions.optionsMap["system.box"].as<std::vector<float>>()[2];
            mobile = programOptions.optionsMap["system.mobile"].as<std::size_t>();
        }
        else 
            vesCRITICAL("UNKNOWN ERROR")


        if(programOptions.optionsMap.count("system.temperature"))
            temperature = programOptions.optionsMap["system.temperature"].as<float>();
        else 
            vesCRITICAL("system.temperature not defined")

        
        if(programOptions.optionsMap.count("system.timestep"))
            dt = programOptions.optionsMap["system.timestep"].as<float>();
        else 
            vesCRITICAL("system.timestep not defined")


        if(programOptions.optionsMap.count("system.kappa"))
            kappa = programOptions.optionsMap["system.kappa"].as<float>();
        else 
            vesCRITICAL("system.kappa not defined")


        if(programOptions.optionsMap.count("system.gamma"))
            gamma = enhance::deg_to_rad(programOptions.optionsMap["system.gamma"].as<float>());
        else 
            vesCRITICAL("system.gamma not defined")   
    }

    // output configuration
    {
        if(programOptions.optionsMap.count("output.traj"))
            traj = programOptions.optionsMap["output.traj"].as<std::string>();
        else 
            vesWARNING("output.traj not defined, no trajectory output")   


        if(programOptions.optionsMap.count("output.skip"))
            traj_skip = programOptions.optionsMap["output.skip"].as<std::size_t>();
        else 
        {
            vesWARNING("output.skip not defined, will print every step")
            traj_skip = 1;
        }
    }

    {
        vesLOG("VALUE OVERVIEW")
        vesLOG("general.algorithm   " << algorithm )
        vesLOG("general.interaction " << interaction )
        vesLOG("general.thermostat  " << thermostat )
        vesLOG("system.mobile       " << mobile )
        vesLOG("system.box          " << x << " " << y << " " << z )
        vesLOG("system.timestep     " << dt )
        vesLOG("system.termperature " << temperature )
        vesLOG("system.kappa        " << kappa )
        vesLOG("system.gamma        " << gamma )
        vesLOG("output.traj         " << traj )
        vesLOG("output.skip         " << traj_skip )
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
