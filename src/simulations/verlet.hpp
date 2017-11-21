#pragma once

#include "algorithm.hpp"
#include <iostream>
#include <tbb/tbb.h>



class Verlet
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