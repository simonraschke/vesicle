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
        float pre_sum = 0;
        for(std::size_t j = 0; j<i; ++j)
            pre_sum += algorithm->getInteraction().value(particles[i],particles[j]);

        auto current = energy_sum.load();
        while (!energy_sum.compare_exchange_weak(current, current + pre_sum))
            current = energy_sum.load();
        
    });
    return !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
}


float System::kineticEnergy() const
{   
    float value = PARALLEL_ACCUMULATE(particles,[&](float f, auto& p){ assert(p); return f + 0.5f*p->velocity().squaredNorm(); });
    return !std::isnan(value) ? value : throw std::runtime_error("kinetic Energy is NAN");
}



void System::addTime(float t)
{
    time_elapsed += t;
}



float System::getTime() const
{
    return time_elapsed;
}