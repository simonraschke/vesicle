#include "system.hpp"



void System::clear()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    algorithm.reset(nullptr);
    thermostat.reset(nullptr);
    trajectory_writer.reset(nullptr);
    particles.clear();
    assert(!algorithm);
    assert(!thermostat);
    assert(!trajectory_writer);
}



const PARTICLERANGE& System::getParticles() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return particles;
}



const std::unique_ptr<Algorithm>& System::getAlgorithm() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(algorithm);
    return algorithm;
}



const std::unique_ptr<Interaction>& System::getInteraction() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(algorithm);
    assert(algorithm->getInteraction());
    return algorithm->getInteraction();
}



const std::unique_ptr<Thermostat>& System::getThermostat() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return thermostat;
}




const std::unique_ptr<TrajectoryWriter>& System::getTrajectoryWriter() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return trajectory_writer;
}



float System::potentialEnergy() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    std::atomic<float> energy_sum = {0.f};
    tbb::parallel_for(std::size_t(0),particles.size(), std::size_t(1), [&](const std::size_t& i)
    {
        float pre_sum = 0;
        for(std::size_t j = 0; j<i; ++j)
        {
            assert(algorithm);
            assert(algorithm->getInteraction());
            assert(particles[i]);
            assert(particles[j]);
            pre_sum += algorithm->getInteraction()->translation(particles[i],particles[j]) + algorithm->getInteraction()->rotation(particles[i],particles[j]);
        }

        auto current = energy_sum.load();
        while (!energy_sum.compare_exchange_weak(current, current + pre_sum))
            current = energy_sum.load();
        
    });
    return !std::isnan(energy_sum.load()) ? energy_sum.load() : throw std::runtime_error("potential Energy is NAN");
}


float System::kineticEnergy() const
{   
    vesDEBUG(__PRETTY_FUNCTION__)
    float value = PARALLEL_REDUCE(float,particles,[&](float f, auto& p){ assert(p); return f + 0.5f*p->velocity().squaredNorm(); });
    return !std::isnan(value) ? value : throw std::runtime_error("kinetic Energy is NAN");
}



void System::addTime(float t)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    time_elapsed += t;
}



float System::getTime() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return time_elapsed;
}