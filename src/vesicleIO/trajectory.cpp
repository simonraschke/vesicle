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

#include "trajectory.hpp"



TrajectoryRWBase::TrajectoryRWBase()
    : working_dir(boost::filesystem::current_path())
{
    vesDEBUG(__PRETTY_FUNCTION__)
}



TrajectoryRWBase::~TrajectoryRWBase()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(FILE.is_open())
    {
        FILE.close();
    }
}



void TrajectoryRWBase::setPath(PATH path)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << path)
    if(!file_path && FILE.is_open())
    {
        FILE.close();
    }
    
    file_path = std::make_unique<PATH>(boost::filesystem::system_complete(working_dir/path));
    if( !boost::filesystem::exists(*file_path) )
        throw std::runtime_error("path does not exist: " + file_path->generic_string());
    assert(boost::filesystem::exists(*file_path));

    FILE.open(*file_path);
}



TrajectoryRWBase::PATH TrajectoryRWBase::getFilePath() const
{
    return *file_path;
}



TrajectoryRWBase::PATH TrajectoryRWBase::getWorkingDir() const
{
    return working_dir;
}



bool TrajectoryRWBase::isOpen() const
{
    return FILE.is_open();
}



bool TrajectoryRWBase::isEOF() const
{
    return FILE.eof();
}



void TrajectoryRWBase::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}



void TrajectoryRWBase::setAnisotropic(bool b)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << b)
    anisotropic = b;
    // makeStartFileVMD();
}



void TrajectoryWriter::setAnisotropic(bool b)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << b)
    anisotropic = b;
    // makeStartFileVMD();
}



void TrajectoryReader::clearAllFrames()
{
    frames.clear();
}
