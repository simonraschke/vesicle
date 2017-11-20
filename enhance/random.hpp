#pragma once

#include <random>
#include <eigen3/Eigen/Core>



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



    // template<typename EigenDerived, typename T>
    // EigenDerived random(const T& a, const T& b)
    // {   
    //     // this will generate an eigen vector with random components
    //     // depending on the implementation of random for the
    //     // Eigen::MatrixBase::Scalar type
    //     return EigenDerived().unaryExpr([&](const typename EigenDerived::Scalar& val){ return random<decltype(val)>(a,b);});
    // }
}