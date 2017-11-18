#pragma once

#include <eigen3/Eigen/Geometry>
#include <memory>
#include "particles/particle_base.hpp"


enum class PERIODIC : bool { ON=true, OFF=false};

template<PERIODIC P>
class Box
{
public:
    typedef ParticleInterface::real real;
    typedef ParticleInterface::cartesian cartesian;

    void setLengthX(real);
    void setLengthY(real);
    void setLengthZ(real);

    void unsetLengthX();
    void unsetLengthY();
    void unsetLengthZ();

    real getLengthX() const;
    real getLengthY() const;
    real getLengthZ() const;

    real distance(const cartesian&, const cartesian&);
    real distance(const ParticleInterface&, const ParticleInterface&);

    cartesian scaleDown(cartesian);
    cartesian scaleDown(const ParticleInterface&);

protected:
    void check_for_aligned_box_setup();

private:
    std::unique_ptr<real> x {nullptr}; 
    std::unique_ptr<real> y {nullptr};
    std::unique_ptr<real> z {nullptr};

    std::unique_ptr<Eigen::AlignedBox<real,3>> bounding_box {nullptr};
};






template<PERIODIC P>
void Box<P>::setLengthX(real l)
{
    x = std::make_unique<real>(l);
    assert(x);
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthY(real l)
{
    y = std::make_unique<real>(l);
    assert(y);
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::setLengthZ(real l)
{
    z = std::make_unique<real>(l);
    assert(z);
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::unsetLengthX()
{
    x.reset(nullptr);
    bounding_box.reset(nullptr);
    assert(!x);
    assert(!bounding_box);
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::unsetLengthY()
{
    y.reset(nullptr);
    bounding_box.reset(nullptr);
    assert(!y);
    assert(!bounding_box);
    check_for_aligned_box_setup();
}



template<PERIODIC P>
void Box<P>::unsetLengthZ()
{
    z.reset(nullptr);
    bounding_box.reset(nullptr);
    assert(!z);
    assert(!bounding_box);
    check_for_aligned_box_setup();
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthX() const
{
    assert(x);
    return x;
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthY() const
{
    assert(y);
    return y;
}



template<PERIODIC P>
typename Box<P>::real Box<P>::getLengthZ() const
{
    assert(z);
    return z;
}



template<PERIODIC P>
void Box<P>::check_for_aligned_box_setup()
{
    if( x && y && z )
    {
        assert(x);
        assert(y);
        assert(z);
        bounding_box = std::make_unique<Eigen::AlignedBox<real,3>>( cartesian::Zero(), cartesian(*x, *y, *z) );
        assert(bounding_box);
    }
}



template<>
Box<PERIODIC::ON>::real Box<PERIODIC::ON>::distance(const cartesian& c1, const cartesian& c2)
{
    assert(bounding_box);
    cartesian distance_cartesian;
    distance_cartesian = c2-c1;
    distance_cartesian(0) = distance_cartesian(0) - (*x) * std::round(distance_cartesian(0)/(*x));
    distance_cartesian(1) = distance_cartesian(1) - (*y) * std::round(distance_cartesian(1)/(*y));
    distance_cartesian(2) = distance_cartesian(2) - (*z) * std::round(distance_cartesian(2)/(*z));
    return distance_cartesian.norm();
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const cartesian& c1, const cartesian& c2)
{
    assert(bounding_box);
    return (c2-c1).norm();
}



template<>
Box<PERIODIC::ON>::real Box<PERIODIC::ON>::distance(const ParticleInterface& p1, const ParticleInterface& p2)
{
    return distance(p2.coords(), p1.coords());
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const ParticleInterface& p1, const ParticleInterface& p2)
{
    return distance(p1.coords(),p2.coords());
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(cartesian c)
{
    c(0) = c(0) - (*x) * std::round(c(0)/(*x));
    c(1) = c(1) - (*y) * std::round(c(1)/(*y));
    c(2) = c(2) - (*z) * std::round(c(2)/(*z));
    return c;
}



template<PERIODIC P>
typename Box<P>::cartesian Box<P>::scaleDown(const ParticleInterface& p)
{
    return scaleDown(p.coords());
}