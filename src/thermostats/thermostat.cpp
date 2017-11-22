#include "thermostat.hpp"



void Thermostat::setTarget(PARTICLERANGE* range)
{
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}