#pragma once

#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    float dt = 0.01;
    float x = 7;
    float y = 7;
    float z = 7;

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