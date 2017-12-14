/*  
*   Copyright 2017 Simon Raschke
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

#include "grid.hpp"



GridGeometry::GridGeometry()
{

}



GridGeometry::GridGeometry(std::size_t _x, std::size_t _y, std::size_t _z)
    : x(_x)
    , y(_y)
    , z(_z)
{
    generate();
}


void GridGeometry::generate()
{
    points.reserve(x*y*z);

    for(std::size_t i = 0; i<x; ++i)
        for(std::size_t j = 0; j<y; ++j)
            for(std::size_t k = 0; k<z; ++k)    
            {
                points.emplace_back(cartesian(i,j,k));
            }
}



Geometry::cartesian& GridGeometry::matrixView(const std::size_t& i, const std::size_t& j, const std::size_t& k)
{
    return points[ i*x*y + j*y + k ];
}



const Geometry::cartesian& GridGeometry::matrixView(const std::size_t& i, const std::size_t& j, const std::size_t& k) const
{
    return points[ i*x*y + j*y + k ];
}



void GridGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
}



void GridGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
}