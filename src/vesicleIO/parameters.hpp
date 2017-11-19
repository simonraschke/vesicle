#pragma once

#include <memory>
#include <exception>
#include <cassert>




// the actual Parameters class

struct Parameters
{
    float dt = 0.1;
    
};





struct ParameterDependentComponent
{
    void setParameters(Parameters);
    const Parameters& getParameters() const;
    virtual ~ParameterDependentComponent() = default;

protected:
    ParameterDependentComponent() = default;

private:
    std::unique_ptr<Parameters> parameters {nullptr};
};