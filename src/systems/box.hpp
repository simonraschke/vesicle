#pragma once

#include <eigen3/Eigen/Geometry>
#include <memory>
#include "../particles/particle_base.hpp"


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

    real distance(const ParticleInterface&, const ParticleInterface&);
    real distance(const cartesian&, const cartesian&);

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
    return distance(p1.coords(),p2.coords());
}



template<>
Box<PERIODIC::OFF>::real Box<PERIODIC::OFF>::distance(const ParticleInterface& p1, const ParticleInterface& p2)
{
    return distance(p1.coords(),p2.coords());
}