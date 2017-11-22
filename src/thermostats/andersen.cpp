#include "andersen.hpp"



void AndersenThermostat::apply()
{
    auto& target = (*target_range)[enhance::random<std::size_t>(0,target_range->size()-1)];
    auto velocities = target->velocity();
    float mean = 0.f;
    float sigma = std::sqrt(getParameters().temperature/1.f);
    std::normal_distribution<float> gaussian(mean, sigma);
    float new_velocity = gaussian(enhance::RandomEngine.pseudo_engine);
    auto dimension = enhance::random<unsigned int>(0,2);
    velocities(dimension) = new_velocity;
    target->setVelocity(velocities);
}