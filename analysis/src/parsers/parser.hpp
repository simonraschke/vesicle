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

#include "definitions.hpp"
#include "systems/box.hpp"
#include "vesicleIO/parameters.hpp"
#include "enhance/observer_ptr.hpp"



template<typename T>
struct Parser
    : public Box<PERIODIC::ON>
    , public virtual ParameterDependentComponent
{
    typedef T result_type;

    void setTarget(PARTICLERANGE*);
    virtual void parse() = 0;
    virtual ~Parser() = default;

    result_type result {};

protected:
    Parser() = default;

    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
};



template<typename T>
inline void Parser<T>::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}