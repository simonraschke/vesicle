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

#include <cmath>
#include <cstdint>
#include <sstream>
#include <exception>
#include <eigen3/Eigen/Core>
    


namespace enhance
{
    // number of digits
    template <typename T>
    std::uint16_t numDigits(T number)
    {
        std::uint16_t digits = 0;
        if (number < 0) digits = 1; // remove this line if '-' counts as a digit
        std::uint32_t helper = static_cast<int>(number);
        while (helper) 
        {
            helper /= 10;
            digits++;
        }
        return digits;
    }
    

    
    template<typename DERIVED1, typename DERIVED2>
    constexpr float directed_angle(const DERIVED1& v1, const DERIVED2& v2)
    {
        assert(std::isfinite(std::atan2(v2.normalized()(1), v2.normalized()(0)) - std::atan2(v1.normalized()(1), v1.normalized()(0))));
        return std::atan2(v2.normalized()(1), v2.normalized()(0)) - std::atan2(v1.normalized()(1), v1.normalized()(0));
    }
    

    
    template<typename DERIVED1, typename DERIVED2>
    constexpr float absolute_angle(const DERIVED1& v1, const DERIVED2& v2)
    {
        return std::abs(directed_angle(v1,v2));
    }
    

    
    template<typename DERIVED1, typename DERIVED2>
    constexpr float normalized_angle(const DERIVED1& v1, const DERIVED2& v2)
    {
        const float angle = directed_angle(v1,v2);
        return angle < 0.f ? angle+M_PI : angle;
    }
        
    
    
    template<typename T, typename ENABLER = typename std::enable_if<std::is_floating_point<T>::value>::type>
    constexpr T deg_to_rad(const T& __deg) noexcept
    {
        return __deg*M_PI/180;
    }
    
    
    
    template<typename T, typename ENABLER = typename std::enable_if<std::is_floating_point<T>::value>::type>
    constexpr T rad_to_deg(const T& __rad) noexcept
    {
        return __rad/M_PI*180;
    }
    
    
    
    template<typename T>
    constexpr T sphere_volume(const T& __rad) noexcept
    {
        return M_PI*__rad*__rad*__rad*4.f/3.f;
    }
    
    
    
    template<typename T>
    constexpr T sphere_surface(const T& __rad) noexcept
    {
        return M_PI*__rad*__rad*4.f;
    }
    
    
    
    template<typename T>
    constexpr T circle_area(const T& __rad) noexcept
    {
        return M_PI*__rad*__rad;
    }
    
    
    
    template<typename T>
    constexpr T circle_area_to_radius(const T& __area)
    {
        return std::sqrt( __area/M_PI );
    }
    
    
    
    template<typename T>
    constexpr T cone_volume(const T& __rad, const T& __height) noexcept
    {
        return enhance::circle_area(__rad)*__height/((T)3.0);
    }
}