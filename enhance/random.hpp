#pragma once

#include <random>



namespace enhance
{
    static struct RandomEngineInit
    {
        RandomEngineInit();
        int getSeed() const;
        std::mt19937_64 pseudo_engine {};

    private:
        std::random_device true_engine {};
        const int seed;

    }RandomEngine;


    // call these functions
    template<typename T>
    T random(const T& a, const T& b);
}