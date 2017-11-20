#pragma once

#include <eigen3/Eigen/Geometry>
#include <memory>
#include "particles/particle.hpp"
#include "vesicleIO/parameters.hpp"


enum class PERIODIC : bool { ON=true, OFF=false};

template<PERIODIC P>
class Box
    : virtual public ParameterDependentComponent
{
public:
    typedef Particle::real real;
    typedef Particle::cartesian cartesian;

    void setLengthX(real);
    void setLengthY(real);
    void setLengthZ(real);

    real getLengthX() const;
    real getLengthY() const;
    real getLengthZ() const;

    real distance(const cartesian&, const cartesian&) const;
    real distance(const Particle&, const Particle&) const;

    real squared_distance(const cartesian&, const cartesian&) const;
    real squared_distance(const Particle&, const Particle&) const;

    cartesian scaleDown(cartesian) const ;
    cartesian scaleDown(const Particle&) const;

    bool contains(const cartesian&) const;
    bool contains(const Particle&) const;

    virtual ~Box() = default;

protected:
    using ParameterDependentComponent::mutableAccess;

    void check_for_aligned_box_setup();

private:

    std::unique_ptr<Eigen::AlignedBox<real,3>> bounding_box {nullptr};
};




template<PERIODIC P>
void Box<P>::setLengthX(real l)
{
    mutableAccess().x = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthY(real l)
{
    mutableAccess().y = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthZ(real l)
{
    mutableAccess().z = l;
    check_for_aligned_box_setup();
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthX() const
{
    return getParameters().x;
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthY() const
{
    return getParameters().x;
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthZ() const
{
    return getParameters().z;
}



template<PERIODIC P>
void Box<P>::check_for_aligned_box_setup()
{
    bounding_box.reset(nullptr);
    cartesian vec;
    vec(0) = getParameters().x;
    vec(1) = getParameters().y;
    vec(2) = getParameters().z;
    bounding_box = std::make_unique<Eigen::AlignedBox<real,3>>( cartesian::Zero(), vec );
}



// template<>
// Box<PERIODIC::ON>::real Box<PERIODIC::ON>::squared_distance(const cartesian& c1, const cartesian& c2) const 
// {
//     cartesian distance_cartesian;
//     distance_cartesian = c2-c1;
//     distance_cartesian(0) = distance_cartesian(0) - getParameters().x * std::round(distance_cartesian(0)/(getParameters().x));
//     distance_cartesian(1) = distance_cartesian(1) - getParameters().y * std::round(distance_cartesian(1)/(getParameters().y));
//     distance_cartesian(2) = distance_cartesian(2) - getParameters().z * std::round(distance_cartesian(2)/(getParameters().z));
//     return distance_cartesian.squaredNorm();
// }



// template<>
// Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::squared_distance(const cartesian& c1, const cartesian& c2) const;



// template<>
// Box<PERIODIC::ON>::real Box<PERIODIC::ON>::squared_distance(const Particle& p1, const Particle& p2) const;



// template<>
// Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::squared_distance(const Particle& p1, const Particle& p2) const;



// template<>
// Box<PERIODIC::ON>::real Box<PERIODIC::ON>::distance(const cartesian& c1, const cartesian& c2) const;



// template<>
// Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const cartesian& c1, const cartesian& c2) const;



// template<>
// Box<PERIODIC::ON>::real Box<PERIODIC::ON>::distance(const Particle& p1, const Particle& p2) const;



// template<>
// Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const Particle& p1, const Particle& p2) const;



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(cartesian c) const 
{
    c(0) = c(0) - getParameters().x * std::round(c(0)/(getParameters().x));
    c(1) = c(1) - getParameters().y * std::round(c(1)/(getParameters().y));
    c(2) = c(2) - getParameters().z * std::round(c(2)/(getParameters().z));
    return c;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(const Particle& p) const 
{
    return scaleDown(p.coords());
}