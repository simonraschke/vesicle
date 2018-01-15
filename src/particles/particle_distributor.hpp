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

#include "definitions.hpp"
#include "particle.hpp"
#include "vesicleIO/parameters.hpp"
#include "enhance/random.hpp"
#include "systems/box.hpp"
#include "geometries/grid.hpp"
#include <tbb/parallel_for_each.h>
#include <atomic>
#include <iostream>



struct Distributor
    : public Box<PERIODIC::ON>
    , virtual public ParameterDependentComponent
{
    virtual void operator()(PARTICLERANGE*) = 0;
    
    virtual ~Distributor() = default;

protected:
    Distributor() = default;

    bool conflicting_placement(PARTICLERANGE*, PARTICLERANGE::value_type&);
};



struct RandomDistributor
    : public Distributor
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual void operator()(PARTICLERANGE*) override;

protected:
    cartesian randomCoords() const;
    cartesian randomOrientation() const;
};



struct GridDistributor
    : public Distributor
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual void operator()(PARTICLERANGE*) override;

protected:
};