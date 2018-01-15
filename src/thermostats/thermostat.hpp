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