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

#include "interaction.hpp"



class AngularLennardJones
    : public Interaction
{
public:
    virtual void setup() override;
    
    virtual float potential(const Particle&, const Particle&) const override;
    virtual cartesian force(const Particle&, const Particle&) const override;

    virtual bool isAnisotropic() const override;
    
protected:
    // float chi(const Particle&, const Particle&) const; 

    float kappa {0};
    float a {0};
    float b {0};
    float c {0};
    float sigma {0};
    float epsilon {0};
    float cutoff_rez_sq {0};
};