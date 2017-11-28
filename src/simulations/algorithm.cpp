#include "algorithm.hpp"


void Algorithm::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}



const std::unique_ptr<Interaction>& Algorithm::getInteraction() const
{
    assert(interaction);
    return interaction;
}