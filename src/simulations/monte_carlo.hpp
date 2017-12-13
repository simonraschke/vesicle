#pragma once

#include "algorithm.hpp"
#include "acceptance_adapters/metropolis.hpp"
#include <memory>



class MonteCarlo
    : public Algorithm
{
public: 
    virtual void step(const unsigned long& = 1) override;
    
protected:
    virtual void updateCoords() override;
    virtual void updateForces() override;
    virtual void updateVelocities() override;

    float potentialEnergy(const std::unique_ptr<Particle>&) const;

private:  
};