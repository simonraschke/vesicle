#include "potential_energy.hpp"



const std::unique_ptr<Interaction>& PotentialEnergyParser::getInteraction() const
{
    assert(interaction);
    return interaction;
}



void PotentialEnergyParser::parse()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!interaction)
        throw std::logic_error("interaction was not set");

    std::atomic<float> energy_sum = {0.f};
    tbb::parallel_for(std::size_t(0),target_range->size(), std::size_t(1), [&](const std::size_t& i)
    {
        float pre_sum = 0;
        for(std::size_t j = 0; j<i; ++j)
        {
            assert();
            assert(interaction);
            assert(particles[i]);
            assert(particles[j]);
            pre_sum += interaction->potential(target_range->at(i),target_range->at(j));
        }

        // this works for an std::atomic float, but may busy wait to compare_exchange like hell if too many particles
        auto current = energy_sum.load();
        while (!energy_sum.compare_exchange_weak(current, current + pre_sum))
            current = energy_sum.load();
        
    });

    result = !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
}
