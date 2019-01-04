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

#include "trajectory.hpp"
#include <deque>
#include <map>
#include <algorithm>
#include <cctype>
#include <boost/algorithm/string.hpp>
// #include <hdf5.h>



// Gromacs trajectory reader class
// interface provided by base class
class TrajectoryReaderH5
    : public virtual TrajectoryReader
{
public:
    TrajectoryReaderH5();
    ~TrajectoryReaderH5() = default;

    virtual void setPath(PATH) override;

    // read input stream and safe last frame
    virtual void readAllFrames(bool = true) override;
    virtual void readNextFrame(std::regex) override;

    virtual float getTime() const override;

protected:
    std::vector<std::string> getGroupNames() const;


private:
};