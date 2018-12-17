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

#include "particle.hpp"



class ParticleFrame
    : public Particle
{
public:
    virtual void setCoords(const cartesian&) override;
    virtual void setVelocity(const cartesian&) override;
    virtual void setForce(const cartesian&) override;
    virtual void setOrientation(const cartesian&) override;

    void setOriginCoords(const cartesian&);
    void setOriginOrientation(const cartesian&);

    virtual void setOffset(float, float, float) override;

    virtual std::string name() const override;
    virtual PARTICLETYPE getType() const override;

protected:
    bool coordsSetupDone = false;
    bool orientationSetupDone = false;
    real offsetX = 0.3;
    real offsetY = 0.3;
    real offsetZ = 0.3;
    std::unique_ptr<cartesian> originCoords {nullptr};
    std::unique_ptr<cartesian> originOrientation {nullptr};
};