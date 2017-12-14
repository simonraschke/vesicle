#pragma once

#include "acceptance_adapter.hpp"
#include <cmath>



class Metropolis
    : public AcceptanceAdapter
{
public:
    virtual bool isValid(float) const override;
};



inline bool Metropolis::isValid(float energy_difference) const
{
    #ifndef NDEBUG
        const float tem = getParameters().temperature;
        const float exr = std::exp(-energy_difference/getParameters().temperature);
        const float ran = enhance::random<float>(0.0,1.0);
        const bool acc = exr > ran;
        vesLOG("energy difference: "<< energy_difference << "  temp: " << tem << "  exp: " << exr <<"  random: " << ran << "  accepted: " << std::boolalpha << acc )
        assert( energy_difference < 0.f ? acc : true);
        return acc;
    #else
        return energy_difference < 0.f ? true : std::exp(-energy_difference/getParameters().temperature) > enhance::random<float>(0.0,1.0);
    #endif
    
}