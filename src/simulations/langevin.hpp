#pragma once

#include <iostream>
#include "algorithm.hpp"



class Langevin
    : public Algorithm
{
public:
    
    virtual void step(const unsigned long& = 1) override;
    
protected:
    virtual void updateForces() override;
    virtual void updateCoords() override;
    virtual void updateOrientations() override;

private:  
};