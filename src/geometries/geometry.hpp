#pragma once

#include "definitions.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>



struct Geometry
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual ~Geometry() = default;
    virtual void generate() = 0;
    virtual void scale(const cartesian&) = 0;
    virtual void shift(const cartesian&) = 0;

    std::vector<cartesian> points {};

protected:
    Geometry() = default;
};