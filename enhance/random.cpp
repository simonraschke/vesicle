#include "random.hpp"



enhance::RandomEngineInit::RandomEngineInit()
    : seed(true_engine())
{
    pseudo_engine.seed(seed);
}



int enhance::RandomEngineInit::getSeed() const
{
    return seed;
}



template<>
float enhance::random(const float& a, const float& b)
{
    std::uniform_real_distribution<float> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}



template<>
double enhance::random(const double& a, const double& b)
{
    std::uniform_real_distribution<double> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}



template<>
unsigned int enhance::random(const unsigned int& a, const unsigned int& b)
{
    std::uniform_int_distribution<unsigned int> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}