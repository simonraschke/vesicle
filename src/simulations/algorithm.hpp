#pragma once

#include "definitions.hpp"
#include "enhance/observer_ptr.hpp"
#include "vesicleIO/parameters.hpp"
#include "particles/particle.hpp"


class Algorithm
    : public ParameterDependentComponent
{
public:

    // set Parameters
    void setTarget(PARTICLERANGE*);
    // void setParameters(Parameters);

    // execute
    virtual void step(const unsigned long& = 1) = 0;

    // necessary
    virtual ~Algorithm() = default;

protected:
    Algorithm() = default;

    virtual void updateForces() = 0;
    virtual void updateCoords() = 0;
    virtual void updateOrientations() = 0;

    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};

private:  
    // std::unique_ptr<Parameters> parameters {nullptr};

};