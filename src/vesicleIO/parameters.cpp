#include "parameters.hpp"





void ParameterDependentComponent::setParameters(Parameters prms)
{
    parameters = std::make_unique<Parameters>(prms);
}



const Parameters& ParameterDependentComponent::getParameters() const
{
    assert(parameters);
    if(!parameters)
    {
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}