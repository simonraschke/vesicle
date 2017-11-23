#pragma once

#include <vector>
#include <memory>



// DEPRECATED
struct HistoryBuffer
{
    std::unique_ptr<float> time {nullptr};
    std::unique_ptr<float> kineticEnergy {nullptr};
    std::unique_ptr<float> potentialEnergy {nullptr};
};



#include "definitions.hpp"
#include <cmath>
#include <string>
#include <boost/filesystem.hpp>
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
    if(buffer.time)
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