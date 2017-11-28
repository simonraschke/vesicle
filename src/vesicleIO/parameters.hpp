#pragma once

#include "enhance/math_utility.hpp"
#include "program_options.hpp"
#include "definitions.hpp"
#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    // GENERAL
    std::string algorithm {};
    std::string interaction {};
    std::string thermostat {};

    // SYSTEM
    std::size_t mobile {};
    float dt {};
    float x {};
    float y {};
    float z {};
    float temperature {};
    float kappa {};
    float gamma {};

    // OUTPUT
    std::string traj {};
    std::size_t traj_skip {};

    // functions
    void setup();
    ProgramOptions programOptions {};
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