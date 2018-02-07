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
#include "parameters.hpp"
#include "history.hpp"
#include "enhance/observer_ptr.hpp"
#include "enhance/output_utility.hpp"
#include "systems/box.hpp"
#include <boost/filesystem.hpp>
#include <iomanip>
#include <string>


// trajectory reader and writer base class
// provides parameters 
//   working directory,
//   filename and filestream
//   box operations

// TODO make setFilename to setFile und PATH with working dir
struct TrajectoryRWBase
    : public virtual ParameterDependentComponent
    , public virtual Box<PERIODIC::ON>
{
    typedef Particle::cartesian cartesian;
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::fstream FSTREAM;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    // setting path to file and opening filestream
    // will set path relative to working directory
    virtual void setPath(PATH);
    PATH getFilePath() const;
    PATH getWorkingDir() const;
    bool isOpen() const;

    // set particle range to write for Writer derived class
    // TODO might be unnecessary in reader class
    void setTarget(PARTICLERANGE*);

    // are particles and interactions anisotropic
    virtual void setAnisotropic(bool);

    // necessary for base classes
    virtual ~TrajectoryRWBase();
    
protected:
    // provide simulation box functions
    using Box<PERIODIC::ON>::scaleDown;
    using Box<PERIODIC::ON>::getLengthX;
    using Box<PERIODIC::ON>::getLengthY;
    using Box<PERIODIC::ON>::getLengthZ;

    // do not construct directly
    TrajectoryRWBase();

    PATH working_dir;
    std::unique_ptr<PATH> file_path {nullptr};
    FSTREAM FILE {};
    // std::unique_ptr<std::string> filename {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
    bool anisotropic {false};
};



// Writer base class for trajectory
// provides writer specific interface
class TrajectoryWriter
    : public TrajectoryRWBase
{
public:
    virtual ~TrajectoryWriter() = default ;
    virtual void setAnisotropic(bool);

    virtual void write(const HistoryStorage&) = 0;

protected:
    virtual void makeStartFileVMD() const = 0;
    TrajectoryWriter() = default;

    unsigned int skip_counter{1};
};



// Writer base class for trajectory
// provides reader specific interface
class TrajectoryReader
    : public TrajectoryRWBase
{
public:
    typedef std::pair<std::size_t,std::deque<std::string>> Frame;
    // typedef std::map<Frame::first_type,Frame::second_type> FrameMap;
    typedef std::vector<Frame> FrameMap;
    
    virtual ~TrajectoryReader() = default;

    virtual void readAllFrames(bool) = 0;
    virtual void readNextFrame(std::regex) = 0;
    void clearAllFrames();

protected:
    TrajectoryReader() = default;

    void readHeaderLine();
    void readParticleLine();
    void readBottomLine();
    
    FrameMap frames {};

    Frame::first_type frame_counter {0};
    Frame::second_type frame_buffer {};

    std::streampos FILE_pos {0};

private:

};