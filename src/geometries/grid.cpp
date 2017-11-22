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