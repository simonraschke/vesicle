#include "algorithm.hpp"


void Algorithm::setTarget(PARTICLERANGE* range)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    target_range = enhance::make_observer<PARTICLERANGE>(range);

    // energy_matrix_old->resize(target_range->size(),target_range->size());
    // energy_matrix_work->resize(target_range->size(),target_range->size());

    // energy_matrix_old->reserve(target_range->size()*20);
    // energy_matrix_work->reserve(target_range->size()*20);
}



const std::unique_ptr<Interaction>& Algorithm::getInteraction() const
{
    assert(interaction);
    return interaction;
}



const std::unique_ptr<AcceptanceAdapter>& Algorithm::getAcceptance() const
{
    assert(acceptance);
    return acceptance;
}