#include "system.hpp"



void System::clear()
{
    algorithm.reset(nullptr);
    particles.clear();
}



Algorithm& System::getAlgorithm() const
{
    return *algorithm;
}



Thermostat& System::getThermostat() const
{
    return *thermostat;
}




TrajectoryWriter& System::getTrajectoryWriter() const
{
    return *trajectory_writer;
}



float System::potentialEnergy() const
{
    std::atomic<float> energy_sum;
    tbb::parallel_for(std::size_t(0),particles.size(), std::size_t(1), [&](const std::size_t& i)
    {
        for(std::size_t j = 0; j<i; ++j)
        {
            const auto& target1 = particles[i];
            const auto& target2 = particles[j];
            auto value = algorithm->getInteraction().value(target1,target2);
            auto current = energy_sum.load();
            while (!energy_sum.compare_exchange_weak(current, current + value))
            {
                current = energy_sum.load();
            }
        }
    });
    return energy_sum;
}


float System::kineticEnergy() const
{
    return PARALLEL_ACCUMULATE(particles,[&](float f, auto& p){ return f + 0.5f*p->velocity().squaredNorm(); });
}