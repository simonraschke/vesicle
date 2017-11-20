#include "box.hpp"



template<>
Box<PERIODIC::ON>::real Box<PERIODIC::ON>::squared_distance(const cartesian& c1, const cartesian& c2) const 
{
    cartesian distance_cartesian;
    distance_cartesian = c2-c1;
    distance_cartesian(0) = distance_cartesian(0) - getParameters().x * std::round(distance_cartesian(0)/(getParameters().x));
    distance_cartesian(1) = distance_cartesian(1) - getParameters().y * std::round(distance_cartesian(1)/(getParameters().y));
    distance_cartesian(2) = distance_cartesian(2) - getParameters().z * std::round(distance_cartesian(2)/(getParameters().z));
    return distance_cartesian.squaredNorm();
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::squared_distance(const cartesian& c1, const cartesian& c2) const 
{
    return (c2-c1).squaredNorm();
}



template<>
Box<PERIODIC::ON>::real Box<PERIODIC::ON>::squared_distance(const Particle& p1, const Particle& p2) const 
{
    return squared_distance(p2.coords(), p1.coords());
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::squared_distance(const Particle& p1, const Particle& p2) const 
{
    return squared_distance(p1.coords(),p2.coords());
}



template<>
Box<PERIODIC::ON>::real Box<PERIODIC::ON>::distance(const cartesian& c1, const cartesian& c2) const 
{
    return std::sqrt(squared_distance(c1,c2));
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const cartesian& c1, const cartesian& c2) const 
{
    return (c2-c1).norm();
}



template<>
Box<PERIODIC::ON>::real Box<PERIODIC::ON>::distance(const Particle& p1, const Particle& p2) const 
{
    return distance(p2.coords(), p1.coords());
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const Particle& p1, const Particle& p2) const 
{
    return distance(p1.coords(),p2.coords());
}