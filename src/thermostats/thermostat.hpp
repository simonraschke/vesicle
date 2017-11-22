#pragma once

#include "systems/box.hpp"
#include "vesicleIO/parameters.hpp"
#include "definitions.hpp"
#include "enhance/observer_ptr.hpp"



struct Thermostat
    : public Box<PERIODIC::ON>
    , public virtual ParameterDependentComponent
{
    void setTarget(PARTICLERANGE*);

    virtual void apply() = 0;

    virtual ~Thermostat() = default;

protected:
    Thermostat() = default;

    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
};