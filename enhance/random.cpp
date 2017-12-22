#include "random.hpp"



enhance::RandomEngineInit::RandomEngineInit()
    : seed(true_engine())
{
    pseudo_engine.seed(seed);
}



enhance::RandomEngineInit::RandomEngineInit(const int __seed)
    : seed(__seed)
{
    pseudo_engine.seed(seed);
}



int enhance::RandomEngineInit::getSeed() const
{
    return seed;
}



/*
 * specializations of enhance::random<T>(const T,const T)
 */


template<>
float enhance::random(const float a, const float b)
{
    std::uniform_real_distribution<float> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}



template<>
double enhance::random(const double a, const double b)
{
    std::uniform_real_distribution<double> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}



template<>
unsigned int enhance::random(const unsigned int a, const unsigned int b)
{
    std::uniform_int_distribution<unsigned int> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}



template<>
std::size_t enhance::random(const std::size_t a, const std::size_t b)
{
    std::uniform_int_distribution<std::size_t> dist(a,b);
    return dist(enhance::RandomEngine.pseudo_engine);
}