#pragma once

#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    float dt = 0.001;
    float x = 1000;
    float y = 1000;
    float z = 1000;

// protected:
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