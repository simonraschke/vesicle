#pragma once

#include "definitions.hpp"
#include "vesicleIO/parameters.hpp"
#include "enhance/random.hpp"



class AcceptanceAdapter
    : public ParameterDependentComponent
{
public:
    virtual bool isValid(float) const = 0;
    virtual ~AcceptanceAdapter() = default;

protected:
    AcceptanceAdapter() = default;
};