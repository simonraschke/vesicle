#pragma once

#include "enhance/math_utility.hpp"
#include "definitions.hpp"
#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    float dt = 0.001;
    float x = 10;
    float y = 10;
    float z = 10;
    float temperature = 0.01;
    unsigned int trajectory_skip = 10;

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