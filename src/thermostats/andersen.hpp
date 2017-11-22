#pragma once

#include "thermostat.hpp"
#include "enhance/random.hpp"



struct AndersenThermostat
    : public Thermostat
{
    virtual void apply() override;

protected:
};