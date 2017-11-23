#pragma once

#include <memory>
#include <string>
#include <eigen3/Eigen/Core>
#include <tbb/spin_mutex.h>



class Particle
{
public:
    virtual ~Particle() {};
    
    typedef float real;
    typedef Eigen::Matrix<real,3,1,0,3,1> cartesian;

    virtual void setCoords(const cartesian&) = 0;
    virtual void setVelocity(const cartesian&) = 0;
    virtual void setForce(const cartesian&) = 0;
    virtual void setOrientation(const cartesian&) = 0;
    virtual void setCircularVelocity(const cartesian&) = 0;
    virtual void setTorque(const cartesian&) = 0;

    virtual std::string name() const = 0;

    void save();

    void clearCoords();
    void clearVelocity();
    void clearForce();
    void clearOrientation();
    void clearCircularVelocity();
    void clearTorque();
    // void addCoords(const cartesian&);
    void addVelocity(const cartesian&);
    void addForce(const cartesian&);
    void addCircularVelocity(const cartesian&);
    void addTorque(const cartesian&);
    // void addOrientation(const cartesian&);
    const cartesian& coords() const;
    const cartesian& coordsOld() const;
    const cartesian& force() const;
    const cartesian& forceOld() const;
    const cartesian& velocity() const;
    const cartesian& velocityOld() const;
    const cartesian& orientation() const;
    const cartesian& orientationOld() const;
    const cartesian& circularVelocity() const;
    const cartesian& circularVelocityOld() const;
    const cartesian& torque() const;
    const cartesian& torqueOld() const;

    void setMass(float);
    float getMass() const;

protected:
    Particle() = default;

    float mass = {1.0};

    std::unique_ptr<cartesian> currentCoords {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldCoords {std::make_unique<cartesian>()};
    
    std::unique_ptr<cartesian> currentForce {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldForce {std::make_unique<cartesian>()};
    
    std::unique_ptr<cartesian> currentVelocity {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldVelocity {std::make_unique<cartesian>()};

    std::unique_ptr<cartesian> currentOrientation {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldOrientation {std::make_unique<cartesian>()};

    std::unique_ptr<cartesian> currentCircularVelocity {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldCircularVelocity {std::make_unique<cartesian>()};

    std::unique_ptr<cartesian> currentTorque {std::make_unique<cartesian>()};
    std::unique_ptr<cartesian> oldTorque {std::make_unique<cartesian>()};

    tbb::spin_mutex mutex {};
private:
};