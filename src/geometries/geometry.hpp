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
#include "enhance/math_utility.hpp"
#include <cmath>
#include <tbb/tbb.h>
#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Core>
#endif
#if __has_include(<Eigen/Core>)
#include <Eigen/Geometry>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Geometry>
#endif



struct Geometry
{
    typedef PARTICLERANGE::value_type::element_type::cartesian cartesian;

    virtual ~Geometry() = default;
    virtual void generate() = 0;
    virtual void scale(const cartesian&) = 0;
    virtual void shift(const cartesian&) = 0;

    tbb::concurrent_vector<cartesian> points {};

protected:
    Geometry() = default;
};