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

#include "sphere.hpp"



SphereGeometry::SphereGeometry()
{

}



SphereGeometry::SphereGeometry(cartesian _origin, float _radius, std::size_t _size)
    : origin(_origin)
    , radius(_radius)
    , size(_size)
{
    generate();
}


void SphereGeometry::generate()
{
    points.reserve(size);

    if(size == 1)
    {
        points[0] = origin + cartesian::UnitZ();
    }
    else if(size == 2)
    {
        points[0] = origin + cartesian::UnitZ();
        points[1] = origin - cartesian::UnitZ();
    }
    else if(size == 3)
    {
        const Eigen::AngleAxisf rotation (M_PI * 2.f / 3.f, cartesian::UnitY());
        points[0] = origin + cartesian::UnitZ();
        points[1] = rotation * points[0];
        points[2] = rotation * points[1];
    }
    else if(size == 4)
    {
        const Eigen::AngleAxisf rotation_down (enhance::deg_to_rad(109.471221) , cartesian::UnitY());
        const Eigen::AngleAxisf rotation_plane (M_PI * 2.f / 3.f , cartesian::UnitZ());
        points[0] = origin + cartesian::UnitZ();
        points[1] = rotation_down * points[0];
        points[2] = rotation_plane * points[1];
        points[3] = rotation_plane * points[2];
    }
    else if(size == 5)
    {
        const Eigen::AngleAxisf rotation_down (M_PI_2 , cartesian::UnitY());
        const Eigen::AngleAxisf rotation_plane (M_PI * 2.f / 3.f , cartesian::UnitZ());
        points[0] = origin + cartesian::UnitZ();
        points[1] = rotation_down * points[0];
        points[2] = rotation_plane * points[1];
        points[3] = rotation_plane * points[2];
        points[4] = rotation_down * points[3];
    }
    else if(size == 6)
    {
        const Eigen::AngleAxisf rotation_down (M_PI_2 , cartesian::UnitY());
        const Eigen::AngleAxisf rotation_plane (M_PI_2 , cartesian::UnitZ());
        points[0] = origin + cartesian::UnitZ();
        points[1] = rotation_down * points[0];
        points[2] = rotation_plane * points[1];
        points[3] = rotation_plane * points[2];
        points[4] = rotation_plane * points[3];
        points[5] = rotation_down * points[4];
    }
    else if(size == 7)
    {
        const Eigen::AngleAxisf rotation_down (M_PI_2 , cartesian::UnitY());
        const Eigen::AngleAxisf rotation_plane (M_PI * 2.f / 5.f , cartesian::UnitZ());
        points[0] = origin + cartesian::UnitZ();
        points[1] = rotation_down * points[0];
        points[2] = rotation_plane * points[1];
        points[3] = rotation_plane * points[2];
        points[4] = rotation_plane * points[3];
        points[5] = rotation_plane * points[4];
        points[6] = rotation_down * points[5];
    }
    else if(size == 8)
    {
        const Eigen::AngleAxisf rotation_down_initial (M_PI_2 , cartesian::UnitY());
        const Eigen::AngleAxisf rotation_down (M_PI_2 , cartesian::UnitY());
        const Eigen::AngleAxisf rotation_plane (M_PI_2 , cartesian::UnitZ());
        points[0] = rotation_down_initial * (origin + cartesian::UnitZ());
        points[1] = rotation_plane * points[0];
        points[2] = rotation_plane * points[1];
        points[3] = rotation_plane * points[2];
        points[4] = rotation_down * points[3];
        points[5] = rotation_plane * points[4];
        points[6] = rotation_plane * points[5];
        points[7] = rotation_plane * points[6];
    }
    else if(size > 8)
    {
        throw std::logic_error("Sphere of size > 8 not yet implemented");
    }

    // scale to radius
    for(auto& point : points)
        point = point.normalized()*radius;
        
}



void SphereGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
}



void SphereGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
}