#pragma once

#include "vesicleIO/parameters.hpp"
#include "particles/particle_base.hpp"


class Algorithm
{
public:
    typedef std::unique_ptr<ParticleInterface> target_type;

    // set Parameters
    void setTarget(std::initializer_list<target_type>);
    void setParameters(Parameters);

    // execute
    virtual void step(const unsigned long& = 1) = 0;

    // necessary
    virtual ~Algorithm() = default;

protected:
    Algorithm();

    virtual void updateForces() = 0;
    virtual void updateCoords() = 0;
    virtual void updateOrientations() = 0;

    std::initializer_list<target_type> target_range ;
private:  
    std::unique_ptr<Parameters> parameters {nullptr};

};