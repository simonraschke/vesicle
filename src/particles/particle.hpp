/*  
*   Copyright 2017 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#pragma once


#include "definitions.hpp"
#include <memory>
#include <string>
#include <eigen3/Eigen/Core>
#include <tbb/spin_mutex.h>



class ParticleIDGenerator
{
public:
    std::size_t operator()() const
    {
        static std::size_t i = 0;
        // ++i;
        return i++;
    }
// private:
    // std::size_t ID = 0;
};




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

    virtual std::string name() const = 0;

    bool operator==(const Particle &);

    void save();
    void clearCoords();
    void clearVelocity();
    void clearForce();
    void clearOrientation();
    // void addCoords(const cartesian&);
    void addVelocity(const cartesian&);
    void addForce(const cartesian&);
    // void addOrientation(const cartesian&);
    const cartesian& coords() const;
    const cartesian& coordsOld() const;
    const cartesian& force() const;
    const cartesian& forceOld() const;
    const cartesian& velocity() const;
    const cartesian& velocityOld() const;
    const cartesian& orientation() const;
    const cartesian& orientationOld() const;

    void setMass(float);
    float getMass() const;

    const unsigned int ID = ParticleIDGenerator()();

protected:
    Particle() = default;

    float mass = {1.0};

    std::unique_ptr<cartesian> currentCoords {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldCoords {std::make_unique<cartesian>(cartesian::Zero())};
    
    std::unique_ptr<cartesian> currentForce {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldForce {std::make_unique<cartesian>(cartesian::Zero())};
    
    std::unique_ptr<cartesian> currentVelocity {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldVelocity {std::make_unique<cartesian>(cartesian::Zero())};

    std::unique_ptr<cartesian> currentOrientation {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldOrientation {std::make_unique<cartesian>(cartesian::Zero())};


    tbb::spin_mutex mutex {};
private:
};