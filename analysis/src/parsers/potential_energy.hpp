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

#include "parser.hpp"
#include "interactions/angular_lennard_jones.hpp"
#include "interactions/lennard_jones.hpp"
#include <atomic>
#include <tbb/atomic.h>
#include <tbb/parallel_for.h>



class PotentialEnergyParser
    : public Parser<float>
{
public:
    virtual void parse() override;

    template<typename I>
    void setInteraction();
    const std::unique_ptr<Interaction>& getInteraction() const;

private:
    std::unique_ptr<Interaction> interaction {nullptr};
};





template<typename I>
void PotentialEnergyParser::setInteraction()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    interaction = std::make_unique<I>();
    interaction->setParameters(getParameters());
    interaction->setup();
}