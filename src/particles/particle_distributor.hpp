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
};



struct GridDistributor
    : public Distributor
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual void operator()(PARTICLERANGE*) override;

protected:
};