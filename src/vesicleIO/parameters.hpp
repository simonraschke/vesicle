#pragma once

#include "enhance/math_utility.hpp"
#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    float dt = 0.01;
    float x = 12;
    float y = 12;
    float z = 12;
    float temperature = 0.5;
    unsigned int trajectory_skip = 100;

    // AngularLennardJones
    float kappa = 1.0;
    float gamma = enhance::deg_to_rad(11.5);
};





struct ParameterDependentComponent
{
    void setParameters(Parameters);
    const Parameters& getParameters() const;
    virtual ~ParameterDependentComponent() = default;

protected:
    Parameters& mutableAccess();
    ParameterDependentComponent() = default;

private:
    std::unique_ptr<Parameters> parameters {nullptr};
};