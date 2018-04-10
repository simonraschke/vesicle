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

#include "geometry.hpp"



struct SphereGeometry
    :public Geometry
{    
    SphereGeometry();
    SphereGeometry(cartesian, float, std::size_t);

    virtual void generate() override;
    virtual void scale(const cartesian&) override;
    virtual void shift(const cartesian&) override;

    cartesian origin {cartesian::Zero()};
    float radius {1.0};
    std::size_t size {100};

    friend std::ostream& operator<<(std::ostream&, const SphereGeometry&);
};