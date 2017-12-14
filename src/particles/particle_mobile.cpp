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

#include "particle_mobile.hpp"



void ParticleMobile::setCoords(const cartesian& newCoords)
{
    assert(!newCoords.hasNaN());
    assert(currentCoords);
    *currentCoords = newCoords;
}



void ParticleMobile::setVelocity(const cartesian& newVelocity)
{
    assert(!newVelocity.hasNaN());
    assert(currentVelocity);
    *currentVelocity = newVelocity;
}



void ParticleMobile::setForce(const cartesian& newForce)
{
    assert(!newForce.hasNaN());
    assert(currentForce);
    *currentForce = newForce;
}



void ParticleMobile::setOrientation(const cartesian& newOrientation)
{
    assert(!newOrientation.hasNaN());
    assert(currentOrientation);
    *currentOrientation = newOrientation.normalized();
}



std::string ParticleMobile::name() const
{
    return "MOBIL";
}