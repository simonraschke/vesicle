/*  
*   Copyright 2017 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

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



    // template<typename EigenDerived, typename T>
    // EigenDerived random(const T& a, const T& b)
    // {   
    //     // this will generate an eigen vector with random components
    //     // depending on the implementation of random for the
    //     // Eigen::MatrixBase::Scalar type
    //     return EigenDerived().unaryExpr([&](const typename EigenDerived::Scalar& val){ return random<decltype(val)>(a,b);});
    // }
}