#pragma once

#include <vector>
#include <memory>



struct HistoryBuffer
{
    std::unique_ptr<float> time {nullptr};
    std::unique_ptr<float> kineticEnergy {nullptr};
    std::unique_ptr<float> potentialEnergy {nullptr};
};


#include <cmath>

struct HistoryStorage
{
    void flush(HistoryBuffer&);
protected:
    std::vector<float> time {};
    std::vector<float> kineticEnergy {};
    std::vector<float> potentialEnergy {};
};



inline void HistoryStorage::flush(HistoryBuffer& buffer)
{
    if(buffer.time)
    {
        time.emplace_back(*buffer.time);
        
        kineticEnergy.emplace_back(buffer.kineticEnergy ? *buffer.kineticEnergy : NAN);
        potentialEnergy.emplace_back(buffer.potentialEnergy ? *buffer.potentialEnergy : NAN);
    }
    else
        return;
}