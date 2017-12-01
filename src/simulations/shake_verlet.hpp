#pragma once

#include "algorithm.hpp"



class ShakeVerlet
    : public Algorithm
{
public:
    virtual void step(const unsigned long& = 1) override;
    
protected:
    virtual void updateCoords() override;
    virtual void updateForces() override;
    virtual void updateVelocities() override;

    void shake();

private:  

    std::unique_ptr<Eigen::MatrixXf> jacobian {nullptr}; 
    std::unique_ptr<Eigen::MatrixXf> jacobianOld {nullptr}; 
    std::unique_ptr<Eigen::MatrixXf> lagrangian {nullptr};
    std::unique_ptr<Eigen::MatrixXf> lagrangianOld {nullptr};
};