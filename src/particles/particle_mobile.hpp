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

#include "particle.hpp"



class ParticleMobile 
    : public Particle
{
public:
    virtual void setCoords(const cartesian&) override;
    virtual void setVelocity(const cartesian&) override;
    virtual void setForce(const cartesian&) override;
    virtual void setOrientation(const cartesian&) override;

    virtual std::string name() const override;
};