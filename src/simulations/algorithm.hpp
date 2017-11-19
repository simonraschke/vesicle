#pragma once

// #include "enhance/observer_ptr.hpp"
#include "particles/particle_base.hpp"


class Algorithm
{
public:
    typedef std::unique_ptr<ParticleInterface> target_type;

    // set Parameters
    void setTarget(std::initializer_list<target_type>);

    // execute
    virtual void step(const unsigned long& = 1) = 0;

    // necessary
    virtual ~Algorithm() = default;

protected:
    Algorithm() : target_range({}) {};

    std::initializer_list<target_type> target_range ;

private:  

};



void Algorithm::setTarget(std::initializer_list<target_type> list)
{
    target_range = list;
}