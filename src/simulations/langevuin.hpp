#pragma once

#include <iostream>
#include "algorithm.hpp"



class Langevuin
    : public Algorithm
{
public:
    
    virtual void step(const unsigned long& = 1) override;
    
protected:
private:  
};



void Langevuin::step(const unsigned long& steps)
{
    for(unsigned long step = 0; step < steps; ++step)
    {
        for(auto& target : target_range)
        {
            std::cout << target->coords() << std::endl << std::endl;
        }
    }
}