#include "parameters.hpp"





void ParameterDependentComponent::setParameters(Parameters prms)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    parameters = std::make_unique<Parameters>(prms);
}



const Parameters& ParameterDependentComponent::getParameters() const
{
    if(!parameters)
    {
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}



Parameters& ParameterDependentComponent::mutableAccess()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!parameters)
    {
        throw std::invalid_argument("parameters is nullptr");
    }
    return *parameters;
}



// Parameters ParameterDependentComponent::getParameters()
// {
//     assert(parameters);
//     if(!parameters)
//     {
//         throw std::invalid_argument("parameters is nullptr");
//     }
//     return *parameters;
// }