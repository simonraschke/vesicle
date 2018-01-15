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

#define USE_MATH_DEFINES

#include "vesicleIO/gro_reader.hpp"
#include "systems/controller.hpp"
#include <csignal>



int main(int argc, const char *argv[])
{
    // register important signals in Controller base class static member
    // allowing civilized shutdown
    std::signal( SIGHUP,  Controller::signal );
    std::signal( SIGINT,  Controller::signal );
    std::signal( SIGQUIT, Controller::signal );
    std::signal( SIGILL,  Controller::signal );
    std::signal( SIGTRAP, Controller::signal );
    std::signal( SIGABRT, Controller::signal );
    std::signal( SIGIOT,  Controller::signal );
    std::signal( SIGBUS,  Controller::signal );
    std::signal( SIGFPE,  Controller::signal );
    std::signal( SIGKILL, Controller::signal );

    TrajectoryReaderGro reader;
    {
        Parameters prms;
        prms.programOptions.read(argc,argv);
        prms.setup();
        reader.setParameters(prms);
    }
    reader.setPath(reader.getParameters().in_traj_path);
    reader.readAllFrames();

    // for( auto frame : reader.getFrames())
    for( auto line : reader.getFrame(-0).second )
    {
        vesLOG(line);
    }

    for( auto frame : reader.getMatches(reader.getParameters().in_frames) )
    {
        for( auto line : frame.second )
        {
            vesLOG(line);
        }
    }

    return EXIT_SUCCESS;
}