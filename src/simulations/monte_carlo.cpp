#include "monte_carlo.hpp"



void MonteCarlo::step(const unsigned long& steps)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    for(unsigned long step = 0; step < steps; ++step)
    {
        assert(target_range);
        tbb::parallel_for_each(target_range->begin(), target_range->end(), [&](auto& target) 
        {
            assert(target);
            target->save();
        });

        updateCoords();
        updateOrientations();
    }
}



void MonteCarlo::updateCoords()
{
    vesDEBUG(__PRETTY_FUNCTION__)

    assert(target_range);
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

        if( !acceptance->isValid(energy_after-energy_before))
        {
            target->setCoords(target->coordsOld());
        }
    });
}



void MonteCarlo::updateOrientations()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    
    assert(target_range);
    std::for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    {
        assert(target);
        const float da = getParameters().stepwidth_orientation;
        const auto randomVec = Particle::cartesian
        (
            enhance::random<float>(0.f,da),
            enhance::random<float>(0.f,da),
            enhance::random<float>(0.f,da)
        );
        const Eigen::AngleAxisf rotate (da, randomVec);

        const float energy_before = potentialEnergy(target);
        target->setOrientation(rotate * target->orientation());
        const float energy_after = potentialEnergy(target);

        if( !acceptance->isValid(energy_after-energy_before))
        {
            target->setOrientation(target->orientationOld());
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
    std::atomic<float> energy_sum {0.f};
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

    // float energy_sum
    // std::for_each(target_range->begin(), target_range->end(), [&](auto& target) 
    // {
    //     assert(p1);
    //     assert(target);
    //     if(p1==target) return;
    //     assert(getInteraction());
    //     const float energy = getInteraction()->translation(p1,target);

    //     auto current = energy_sum.load();
    //     while (!energy_sum.compare_exchange_weak(current, current + energy))
    //         current = energy_sum.load();
    // });
    return !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
}
