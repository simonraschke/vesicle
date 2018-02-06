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



void Parameters::setup()
{
    if(programOptions.optionsMap.count("general.algorithm"))
        algorithm = programOptions.optionsMap["general.algorithm"].as<std::string>();
    else 
    {
        vesWARNING("general.algorithm not defined, choosing verlet")
        algorithm = "verlet";
    }

    if(programOptions.optionsMap.count("general.acceptance"))
        acceptance = programOptions.optionsMap["general.acceptance"].as<std::string>();
    else 
    {
        vesWARNING("general.acceptance not defined, choosing metropolis")
        acceptance = "metropolis";
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
        if(programOptions.optionsMap.count("system.box.x") 
        && programOptions.optionsMap.count("system.box.y") 
        && programOptions.optionsMap.count("system.box.z")) ++counter;

        if(counter != 2)
        {
            vesCRITICAL("define 2 of the 3: system.mobile, system.density, system.box.{_}!")
        }
        else if(programOptions.optionsMap.count("system.mobile") && programOptions.optionsMap.count("system.density"))
        {
            vesLOG("system.mobile and system.density were set. Assuming cubic box")
            mobile = programOptions.optionsMap["system.mobile"].as<std::size_t>();
            x = std::cbrt(static_cast<float>(mobile)/programOptions.optionsMap["system.density"].as<float>());
            y = std::cbrt(static_cast<float>(mobile)/programOptions.optionsMap["system.density"].as<float>());
            z = std::cbrt(static_cast<float>(mobile)/programOptions.optionsMap["system.density"].as<float>());
        }
        else if(programOptions.optionsMap.count("system.box.x") && programOptions.optionsMap.count("system.density"))
        {
            x = programOptions.optionsMap["system.box.x"].as<float>();
            y = programOptions.optionsMap["system.box.y"].as<float>();
            z = programOptions.optionsMap["system.box.z"].as<float>();
            mobile = std::round( programOptions.optionsMap["system.density"].as<float>() * x * y * z);
        }
        else if(programOptions.optionsMap.count("system.box.x") && programOptions.optionsMap.count("system.mobile"))
        {
            x = programOptions.optionsMap["system.box.x"].as<float>();
            y = programOptions.optionsMap["system.box.y"].as<float>();
            z = programOptions.optionsMap["system.box.z"].as<float>();
            mobile = programOptions.optionsMap["system.mobile"].as<std::size_t>();
        }
        else 
            vesCRITICAL("UNKNOWN ERROR")


        if(programOptions.optionsMap.count("system.temperature"))
            temperature = programOptions.optionsMap["system.temperature"].as<float>();
        else 
            vesCRITICAL("system.temperature not defined")

        
        if(programOptions.optionsMap.count("system.timestep"))
        {
            dt = programOptions.optionsMap["system.timestep"].as<float>();
            if(algorithm == std::string("montecarlo"))
                dt = 1;
        }
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


        if(programOptions.optionsMap.count("system.stepwidth_coordinates"))
            stepwidth_coordinates = programOptions.optionsMap["system.stepwidth_coordinates"].as<float>();
        else 
            vesCRITICAL("system.stepwidth_coordinates not defined")   


        if(programOptions.optionsMap.count("system.stepwidth_orientation"))
            stepwidth_orientation = programOptions.optionsMap["system.stepwidth_orientation"].as<float>();
        else 
            vesCRITICAL("system.stepwidth_orientation not defined")   


        if(programOptions.optionsMap.count("system.cell_min_edge"))
            cell_min_edge = programOptions.optionsMap["system.cell_min_edge"].as<float>();
        else 
            vesCRITICAL("system.cell_min_edge not defined")   


        if(programOptions.optionsMap.count("system.max_cells_dim"))
            max_cells_dim = programOptions.optionsMap["system.max_cells_dim"].as<std::size_t>();
        else 
            vesCRITICAL("system.max_cells_dim not defined")   
    }

    // output configuration
    {
        if(programOptions.optionsMap.count("output.traj"))
        {
            std::vector<std::string> possibilities = {"none","gro"};
            if(std::find( std::begin(possibilities), std::end(possibilities), programOptions.optionsMap["output.traj"].as<std::string>() ) != std::end(possibilities) )
                out_traj = programOptions.optionsMap["output.traj"].as<std::string>();
        }
        else 
            vesWARNING("output.out_traj not defined, no trajectory output")   


        if(programOptions.optionsMap.count("output.skip"))
            out_traj_skip = programOptions.optionsMap["output.skip"].as<std::size_t>();
        else 
        {
            vesWARNING("output.skip not defined, will print every step")
            out_traj_skip = 1;
        }
    }

    // input configuration
    {
        if(programOptions.optionsMap.count("input.traj"))
        {
            in_traj = programOptions.optionsMap["input.traj"].as<std::string>();
            if(    in_traj != std::string("none") 
                && in_traj != std::string("gro") )
            {
                vesWARNING("input.traj=" << in_traj << " format not supported, setting to \"none\"")
                in_traj = "none";
            }
        }
        else 
        {
            vesWARNING("input.traj not defined, setting to \"none\"") 
            in_traj = "none";
        }

        if(programOptions.optionsMap.count("input.path"))
        {
            in_traj_path = programOptions.optionsMap["input.path"].as<boost::filesystem::path>();
            if(in_traj != std::string("none") && !boost::filesystem::exists(in_traj_path) && in_traj_path != std::string("none"))
            {
                vesCRITICAL("input.path=" << in_traj_path << " input trajectory not found, setting to \"none\"")
                in_traj_path = "none";
            }
        }
        else 
        {
            vesWARNING("no path to input trajectory defined. setting input.traj to \"none\"") 
            in_traj_path = "none";
        }

        if(programOptions.optionsMap.count("input.frames"))
            in_frames = programOptions.optionsMap["input.frames"].as<std::string>();
        else 
        {
            vesWARNING("no frame defined to read input from. setting to -1 (last frame)") 
            in_frames = "-1";
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
        vesLOG("input.frames                 " << programOptions.optionsMap["input.frames"].as<std::string>() )
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
