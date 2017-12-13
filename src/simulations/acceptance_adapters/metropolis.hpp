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
    return std::exp(-energy_difference/getParameters().temperature) > enhance::random<float>(0.0,1.0);
}