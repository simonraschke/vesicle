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

#include "algorithm.hpp"



void Algorithm::setup()
{
    vesDEBUG(__PRETTY_FUNCTION__ << "nothing to setup")
}



void Algorithm::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);

    // energy_matrix_old->resize(target_range->size(),target_range->size());
    // energy_matrix_work->resize(target_range->size(),target_range->size());

    // energy_matrix_old->reserve(target_range->size()*20);
    // energy_matrix_work->reserve(target_range->size()*20);
}



Algorithm::cell_container_type& Algorithm::getCells()
{
    vesCRITICAL("Algorithm base class cant return cells");
    // return cell_container_type();
}



std::size_t Algorithm::getCurrentStep() const
{
    vesCRITICAL(__PRETTY_FUNCTION__ << " not implemented for Algorithm base class");
    return 0;
}



const std::unique_ptr<Interaction>& Algorithm::getInteraction() const
{
    assert(interaction);
    return interaction;
}



const std::unique_ptr<AcceptanceAdapter>& Algorithm::getAcceptance() const
{
    assert(acceptance);
    return acceptance;
}




void Algorithm::updateCoords()
{
    throw std::logic_error("called virtual function of base class");
}



void Algorithm::updateOrientations()
{
    throw std::logic_error("called virtual function of base class");
}



void Algorithm::updateForces()
{
    throw std::logic_error("called virtual function of base class");
}



void Algorithm::updateVelocities()
{
    throw std::logic_error("called virtual function of base class");
}