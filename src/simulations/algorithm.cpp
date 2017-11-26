#include "algorithm.hpp"


void Algorithm::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}



Interaction& Algorithm::getInteraction() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return *interaction;
}