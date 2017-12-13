#include "monte_carlo.hpp"



void MonteCarlo::step(const unsigned long& steps)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);
        updateCoords();
        updateForces();
        updateVelocities();
    }
}



void MonteCarlo::updateCoords()
{
    std::for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        const float dc = getParameters().stepwidth_coordinates;
        const auto movement = Particle::cartesian
        (
            enhance::random<float>(-dc,dc),
            enhance::random<float>(-dc,dc),
            enhance::random<float>(-dc,dc)
        );
        const float energy_before = potentialEnergy(target);
        target->setCoords(target->coords()+movement);
        const float energy_after = potentialEnergy(target);

        if(acceptance->isValid(energy_after-energy_before))
        {
            
        }
    });
}



void MonteCarlo::updateForces()
{

}



void MonteCarlo::updateVelocities()
{

}



float MonteCarlo::potentialEnergy(const std::unique_ptr<Particle>& p1) const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    std::atomic<float> energy_sum;
    tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(p1);
        assert(target);
        if(p1==target) return;
        assert(getInteraction());
        const float energy = getInteraction()->translation(p1,target);

        auto current = energy_sum.load();
        while (!energy_sum.compare_exchange_weak(current, current + energy))
            current = energy_sum.load();
    });
    return !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
}
