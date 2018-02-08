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

#include "vesicleIO/parameters.hpp"
#include "vesicleIO/trajectory.hpp"
#include <exception>



// translate string to object
// virtual base class
template<typename target_type>
class Translator
    : public virtual ParameterDependentComponent
{
public:
    virtual ~Translator() = default;


protected:
    inline virtual target_type operator()(std::string) const { throw std::logic_error("Translator usage wrong"); };
    inline virtual target_type operator()(std::string,std::string) const { throw std::logic_error("Translator usage wrong"); };
    inline virtual target_type operator()(TrajectoryReader::Frame) { throw std::logic_error("Translator usage wrong"); };
    Translator() = default;

private:
};