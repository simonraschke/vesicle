#include "algorithm.hpp"


void Algorithm::setTarget(PARTICLERANGE* range)
{
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}



Interaction& Algorithm::getInteraction() const
{
    return *interaction;
}