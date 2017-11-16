#pragma once

#include <memory>
#include <eigen3/Eigen/Core>



class ParticleInterface
{
public:
    typedef Eigen::Vector3f cartesian;

    virtual void updateCoords(cartesian&&) = 0;
    virtual void updateOrientation(cartesian&&) = 0;
    const cartesian& coords() const;
    const cartesian& orientation() const;
    const cartesian& coordsOld() const;
    const cartesian& orientationOld() const;

protected:
    ParticleInterface() = default;
    virtual ~ParticleInterface() {};

    std::unique_ptr<cartesian> currentCoords {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldCoords {std::make_unique<cartesian>()};

    std::unique_ptr<cartesian> currentOrientation {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldOrientation {std::make_unique<cartesian>()};
private:
};