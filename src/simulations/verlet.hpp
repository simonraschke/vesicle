#pragma once

#include "algorithm.hpp"



class Verlet
    : public Algorithm
{
public: 
    virtual void step(const unsigned long& = 1) override;
    
protected:
    virtual void updateCoords() override;
    virtual void updateForces() override;
    virtual void updateVelocities() override;

private:  
};