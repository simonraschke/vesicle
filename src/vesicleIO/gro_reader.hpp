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



// Gromacs trajectory reader class
// interface provided by base class
class TrajectoryReaderGro
    : public virtual TrajectoryReader
{
public:
    typedef std::pair<std::size_t,std::deque<std::string>> Frame;
    typedef std::map<Frame::first_type,Frame::second_type> FrameMap;

    TrajectoryReaderGro();
    ~TrajectoryReaderGro() = default;

    // read input stream and safe last frame
    virtual void readAllFrames() override;

    // return last frame 
    // readAllFrames must called beforehand
    Frame getFrame(long long = -1) const;
    FrameMap getMatches(std::regex) const;
    const FrameMap& getFrames() const;

    std::size_t numParticles() const;
    bool isAnisotropic() const;

    static std::map<std::string,std::string> particleLineTokens(std::string line);

protected:
    FrameMap frames {};

private:
};