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
#include <boost/algorithm/string.hpp>



// Gromacs trajectory reader class
// interface provided by base class
class TrajectoryReaderGro
    : public virtual TrajectoryReader
{
public:
    typedef std::pair<std::size_t,std::deque<std::string>> Frame;

    TrajectoryReaderGro();
    ~TrajectoryReaderGro() = default;

    // read input stream and safe last frame
    virtual void readAllFrames() override;

    // return last frame 
    // readAllFrames must called beforehand
    Frame getFrame(long long = -1) const;
    std::map<Frame::first_type,Frame::second_type> getMatches(std::regex) const;
    const std::map<Frame::first_type,Frame::second_type>& getFrames() const;

protected:
    // std::deque<std::string> last_frame {};

    std::map<std::size_t,std::deque<std::string>> frames {};

private:
};