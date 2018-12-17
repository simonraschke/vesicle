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
    if(coordsSetupDone && GLOBAL::getInstance().status == GLOBAL::RUNNING)
    {
        assert(currentCoords);
        if( boundingBox_was_set )
        {
            assert(bounding_box);
            if(bounding_box->contains(newCoords))
            {
                *currentCoords = newCoords;
            }
        }
        else
        {
            if( std::abs(newCoords(0) - (*originCoords)(0)) < offsetX && std::abs(newCoords(1) - (*originCoords)(1)) < offsetY && std::abs(newCoords(2) - (*originCoords)(2)) < offsetZ )
            {
                *currentCoords = newCoords;
            }
        }
        // if( (newCoords-(*originCoords)).squaredNorm() < 0.3f*0.3f )
        // {
        //     *currentCoords = newCoords;
        // }
    }
    else
    {
        originCoords = std::make_unique<cartesian>(newCoords);
        currentCoords = std::make_unique<cartesian>(newCoords);
        coordsSetupDone = true;
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
    if(orientationSetupDone && GLOBAL::getInstance().status == GLOBAL::RUNNING)
    {
        assert(currentOrientation);
        if( std::acos(newOrientation.normalized().dot(originOrientation->normalized())) < M_PI_4/2 )
        {
            *currentOrientation = newOrientation;
        }
        currentOrientation->normalize();
        assert(currentOrientation->squaredNorm()-1.f < 1e-3);
    }
    else
    {
        originOrientation = std::make_unique<cartesian>(newOrientation);
        currentOrientation = std::make_unique<cartesian>(newOrientation);
        assert(originOrientation->squaredNorm()-1.f < 1e-3);
        orientationSetupDone = true;
    }
}



void ParticleFrame::setOriginCoords(const cartesian& origin)
{
    assert(!origin.hasNaN());
    assert(originCoords);
    originCoords = std::make_unique<cartesian>(origin);
    assert(originOrientation->squaredNorm()-1.f < 1e-3);
    coordsSetupDone = true;
}



void ParticleFrame::setOriginOrientation(const cartesian& origin)
{
    assert(!origin.hasNaN());
    assert(originOrientation);
    originOrientation = std::make_unique<cartesian>(origin);
    assert(originOrientation->squaredNorm()-1.f < 1e-3);
    orientationSetupDone = true;
}



void ParticleFrame::setOffset(float x, float y, float z)
{
    offsetX = x;
    offsetY = y;
    offsetZ = z;
}



std::string ParticleFrame::name() const
{
    return "FRAME";
}



PARTICLETYPE ParticleFrame::getType() const
{
    return PARTICLETYPE::FRAME;
}