#pragma once

#include <cmath>
#include <cstdint>


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
        
    
    
    template<typename T>
    constexpr T deg_to_rad(const T& __deg) noexcept
    {
        return __deg*M_PI/180;
    }
    
    
    
    template<typename T>
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