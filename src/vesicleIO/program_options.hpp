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

#include "definitions.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <csignal>
#include <regex>


struct ProgramOptions
{
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    // static void read_from_file(boost::program_options::options_description&, boost::program_options::variables_map&);

    // void read(int, const char* []);

    // boost::program_options::variables_map optionsMap {};
};