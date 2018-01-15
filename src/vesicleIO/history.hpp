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

#include <memory>



// DEPRECATED
struct HistoryBuffer
{
    std::unique_ptr<float> time {nullptr};
    std::unique_ptr<float> kineticEnergy {nullptr};
    std::unique_ptr<float> potentialEnergy {nullptr};
};



// #include "definitions.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <eigen3/Eigen/Core>



// DEPRECATED
struct HistoryStorage
{
    void flush(HistoryBuffer&);
    void dumpToFile(std::string);

    const std::vector<float>& getTime() const;
protected:
    std::vector<float> times {};
    std::vector<float> kineticEnergy {};
    std::vector<float> potentialEnergy {};
};



inline void HistoryStorage::flush(HistoryBuffer& buffer)
{
    if(likely(buffer.time))
    {
        times.emplace_back(*buffer.time);
        
        kineticEnergy.emplace_back(buffer.kineticEnergy ? *buffer.kineticEnergy : NAN);
        potentialEnergy.emplace_back(buffer.potentialEnergy ? *buffer.potentialEnergy : NAN);
    }
    else
        return;
}



inline void HistoryStorage::dumpToFile(std::string n)
{
    assert(times.size() == kineticEnergy.size());
    assert(times.size() == potentialEnergy.size());
    boost::filesystem::ofstream FILE;
    FILE.open(n);

    Eigen::MatrixXf mat(times.size(),3);
    mat.col(0) = Eigen::VectorXf::Map(&times[0],times.size(),1);
    mat.col(1) = Eigen::VectorXf::Map(&kineticEnergy[0],kineticEnergy.size(),1);
    mat.col(2) = Eigen::VectorXf::Map(&potentialEnergy[0],potentialEnergy.size(),1);

    FILE << mat;
    FILE.close();
}



inline const decltype(HistoryStorage::times)& HistoryStorage::getTime() const
{
    return times;
}