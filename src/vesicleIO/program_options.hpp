#pragma once

#include "definitions.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <csignal>


struct ProgramOptions
{
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    static void read_from_file(boost::program_options::options_description&, boost::program_options::variables_map&);

    void read(int, const char* []);

    boost::program_options::variables_map optionsMap {};
};