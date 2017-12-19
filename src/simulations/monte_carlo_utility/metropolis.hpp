/*  
*   Copyright 2017 Simon Raschke
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
        vesDEBUG("energy difference: "<< energy_difference << "  temp: " << tem << "  exp: " << exr <<"  random: " << ran << "  accepted: " << std::boolalpha << acc )
        assert( energy_difference < 0.f ? acc : true);
        return acc;
    #else
        return energy_difference < 0.f ? true : std::exp(-energy_difference/getParameters().temperature) > enhance::random<float>(0.0,1.0);
    #endif
    
}