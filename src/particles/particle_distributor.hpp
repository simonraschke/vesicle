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
#include "vesicleIO/gro_reader.hpp"
#include "enhance/random.hpp"
#include "systems/box.hpp"
#include "geometries/grid.hpp"
#include "geometries/sphere.hpp"
#include "geometries/sphere_grid.hpp"
#include "geometries/plane.hpp"
#include "enhance/incremental_number_gen.hpp"
#include <tbb/parallel_for_each.h>
#include <atomic>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>



struct Distributor
    : public Box<PERIODIC::ON>
    , virtual public ParameterDependentComponent
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual void operator()(PARTICLERANGE*) = 0;
    cartesian randomCoords() const;
    cartesian randomOrientation() const;
    
    virtual ~Distributor() = default;

protected:
    Distributor() = default;

    bool conflicting_placement(PARTICLERANGE*, const PARTICLERANGE::value_type&) const;
};



struct RandomDistributor
    : public Distributor
{
    virtual void operator()(PARTICLERANGE*) override;
    
};



struct GridDistributor
    : public Distributor
{
    virtual void operator()(PARTICLERANGE*) override;

protected:
};



struct OsmoticSystemDistributor
    : public Distributor
{
    virtual void operator()(PARTICLERANGE*) override;

protected:
};



struct TrajectoryDistributorGro
    : public Distributor
{
    virtual void operator()(PARTICLERANGE*) override;

protected:
    typedef std::result_of<decltype(&TrajectoryReaderGro::particleLineTokens)(TrajectoryReaderGro,std::string)>::type tokens_type;
    // typedef decltype(&TrajectoryReaderGro::particleLineTokens) tokens_type;
    virtual void setupAnisotropicParticle(const tokens_type&, const tokens_type&, PARTICLERANGE::value_type::element_type&);
    virtual void setupIsotropicParticle(const tokens_type&,  PARTICLERANGE::value_type::element_type&);
};



struct FrameGuidedGridDistributor
    : public Distributor
{
    virtual void operator()(PARTICLERANGE*) override;
protected:
};



struct FrameGuidedPlaneDistributor
    : public Distributor
{
    virtual void operator()(PARTICLERANGE*) override;
protected:
};