#pragma once

#include "enhance/math_utility.hpp"
#include "definitions.hpp"
#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    float dt = 0.002;
    float x = 5;
    float y = 5;
    float z = 5;
    float temperature = 0.5;
    unsigned int trajectory_skip = 1;

    // AngularLennardJones
    float kappa = 3.0;
    float gamma = enhance::deg_to_rad<float>(20);
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