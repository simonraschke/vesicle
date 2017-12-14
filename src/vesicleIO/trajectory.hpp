/*  
*   Copyright 2017 Simon Raschke
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
#include "systems/box.hpp"
#include <boost/filesystem.hpp>
#include <iomanip>
#include <string>



struct TrajectoryWriter
    : public virtual ParameterDependentComponent
    , public virtual Box<PERIODIC::ON>
{
    typedef Particle::cartesian cartesian;
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    // virtual void read(PATH) = 0;
    virtual void write(const HistoryStorage&) = 0;
    virtual void setFilename(std::string) = 0;

    void setTarget(PARTICLERANGE*);
    void setAnisotropic(bool);

    virtual ~TrajectoryWriter();

protected:
    using Box<PERIODIC::ON>::scaleDown;
    using Box<PERIODIC::ON>::getLengthX;
    using Box<PERIODIC::ON>::getLengthY;
    using Box<PERIODIC::ON>::getLengthZ;

    virtual void makeStartFileVMD() const = 0;

    TrajectoryWriter();

    PATH working_dir;
    OFSTREAM FILE {};
    std::unique_ptr<std::string> filename {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
    unsigned int skip_counter{1};
    bool anisotropic {false};
};



struct TrajectoryWriterGro
    : public TrajectoryWriter
{
    TrajectoryWriterGro();

    virtual void setFilename(std::string) override;
    virtual void write(const HistoryStorage&) override;   

protected:
    virtual void makeStartFileVMD() const override;
};
