#pragma once

#include <memory>
#include <eigen3/Eigen/Core>



class Particle
{
public:
    virtual ~Particle() {};
    
    typedef float real;
    typedef Eigen::Matrix<real,3,1,0,3,1> cartesian;

    virtual void updateCoords(const cartesian&) = 0;
    virtual void updateOrientation(const cartesian&) = 0;

    const cartesian& coords() const;
    const cartesian& orientation() const;
    const cartesian& coordsOld() const;
    const cartesian& orientationOld() const;

protected:
    Particle() = default;

    std::unique_ptr<cartesian> currentCoords {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldCoords {std::make_unique<cartesian>()};

    std::unique_ptr<cartesian> currentOrientation {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldOrientation {std::make_unique<cartesian>()};
private:
};