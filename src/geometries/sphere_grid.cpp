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



#include "sphere_grid.hpp"




SphereGridGeometry::SphereGridGeometry()
{

}



SphereGridGeometry::SphereGridGeometry(std::size_t _edge, float _radius, std::size_t _points)
    : x(_edge)
    , y(_edge)
    , z(_edge)
    , radius(_radius)
    , spherepoints(_points)
{
    generate();
}



SphereGridGeometry::SphereGridGeometry(std::size_t _x, std::size_t _y, std::size_t _z, float _radius, std::size_t _points)
    : x(_x)
    , y(_y)
    , z(_z)
    , radius(_radius)
    , spherepoints(_points)
{
    generate();
}



void SphereGridGeometry::generate()
{
    // make the grid
    points.reserve(x*y*z);
    for(std::size_t i = 0; i<x; ++i)
        for(std::size_t j = 0; j<y; ++j)
            for(std::size_t k = 0; k<z; ++k) 
                points.emplace_back(cartesian(i,j,k));

    set_sphere_points_new();
}



void SphereGridGeometry::set_sphere_points_new()
{
    spheres.clear();
    // set spheres around the grid points
    for(const auto& point : points)
        spheres.emplace_back(point, radius, spherepoints);
}



void SphereGridGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
    set_sphere_points_new();
}



void SphereGridGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
    set_sphere_points_new();
}