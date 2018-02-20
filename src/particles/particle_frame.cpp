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

#include "particle_frame.hpp"



void ParticleFrame::setCoords(const cartesian& newCoords)
{
    assert(!newCoords.hasNaN());
    assert(currentCoords);
    assert(originCoords);
    if( (newCoords-(*originCoords)).norm() < 0.5f )
    {
        *currentCoords = newCoords;
    }
}



void ParticleFrame::setVelocity(const cartesian& newVelocity)
{
    assert(!newVelocity.hasNaN());
    assert(currentVelocity);
    *currentVelocity = newVelocity;
}



void ParticleFrame::setForce(const cartesian& newForce)
{
    assert(!newForce.hasNaN());
    assert(currentForce);
    *currentForce = newForce;
}



void ParticleFrame::setOrientation(const cartesian& newOrientation)
{
    assert(!newOrientation.hasNaN());
    assert(currentOrientation);
    assert(originOrientation);
    if( std::acos(newOrientation.normalized().dot(originOrientation->normalized()) < M_PI_4/2.f) )
    {
        *currentOrientation = newOrientation;
    }
    currentOrientation->normalize();
    assert(currentOrientation->squaredNorm()-1.f < 1e-3);
}



void ParticleFrame::setOriginCoords(const cartesian& origin)
{
    assert(!origin.hasNaN());
    assert(originCoords);
    originCoords = std::make_unique<cartesian>(origin);
    assert(originOrientation->squaredNorm()-1.f < 1e-3);
}



void ParticleFrame::setOriginOrientation(const cartesian& origin)
{
    assert(!origin.hasNaN());
    assert(originOrientation);
    originOrientation = std::make_unique<cartesian>(origin);
    assert(originOrientation->squaredNorm()-1.f < 1e-3);
}



std::string ParticleFrame::name() const
{
    return "FRAME";
}



PARTICLETYPE ParticleFrame::getType() const
{
    return PARTICLETYPE::FRAME;
}