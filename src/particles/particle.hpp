/*  
*   Copyright 2017-2018 Simon Raschke
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
#include "enhance/incremental_number_gen.hpp"
#include <memory>
#include <string>
#include <type_traits>
#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#include <Eigen/Geometry>
#elif __has_include(<eigen3/Eigen/Core>)
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#endif

#include <tbb/spin_mutex.h>



class ParticleIDGenerator : public enhance::IncrementalNumberGenerator<ParticleIDGenerator> {};



enum PARTICLETYPE {MOBILE, FRAME, OSMOTIC, UNDEFINED=-1};


/*
CONCEPT

- all particle classes inherit from Particle
- simulation routine classes are written for containers of type 
  container<std::unique_ptr<Particle>> and must work woth them
  as the overriden interface of derived permits
- members are implemented as pointers since they are smaller than
  Eigen::Vector3 classes without performance loss
- member verctors are initialized zero, mass is one
- set functions are pure virtual and have to be derived

*/


class Particle
{
public:
    // inheritance reuired
    virtual ~Particle() {};

    friend std::ostream& operator<<(std::ostream& os, const Particle&);
    
    // setting the floating point precision
    typedef float real;

    // setting the vector type
    typedef Eigen::Matrix<real,3,1,0,3,1> cartesian;

    // pure virtual setter functions
    // derived overrides must check pointers berfore
    // accessing underlying pointer
    virtual void setCoords(const cartesian&) = 0;
    virtual void setVelocity(const cartesian&) = 0;
    virtual void setForce(const cartesian&) = 0;
    virtual void setOrientation(const cartesian&) = 0;

    // derived must set name for VMD
    virtual std::string name() const = 0;

    // compare by std::addressof comparison
    bool operator==(const Particle &);

    // will set all old member versions to the actual ones
    void save();

    // reset members to 0
    void clearCoords();
    void clearVelocity();
    void clearForce();
    void clearOrientation();

    // adding to velocity and force for verlet purposes
    void addVelocity(const cartesian&);
    void addForce(const cartesian&);

    // return vectors 
    // in DEBUG mode: 
    // - check pointers beforehand
    // - check for NaN's
    const cartesian& coords() const;
    const cartesian& coordsOld() const;
    const cartesian& force() const;
    const cartesian& forceOld() const;
    const cartesian& velocity() const;
    const cartesian& velocityOld() const;
    const cartesian& orientation() const;
    const cartesian& orientationOld() const;

    // setter and getter of mass member
    void setMass(float);
    real getMass() const;

    // virtual function for guiding elements to overwrite constraints
    virtual void setOffset(float, float, float);

    // virtual function for particle types to use
    void setBoundingBox(Eigen::AlignedBox<Particle::real,3>);

    // id will be set automatically
    // id is only necessary to generate useful energy matrices
    // for MonteCarlo simulation routine
    const unsigned int ID = ParticleIDGenerator()();

    // a simple type getter, for easier usage in ParticleSimple
    virtual PARTICLETYPE getType() const = 0;

protected:
    // derive or dont use
    Particle() = default;

    // members
    real mass = {1.0};
    std::unique_ptr<cartesian> currentCoords {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldCoords {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> currentOrientation {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldOrientation {std::make_unique<cartesian>(cartesian::Zero())};

    std::unique_ptr<Eigen::AlignedBox<Particle::real,3>> bounding_box {nullptr};
    bool boundingBox_was_set = false;

    std::unique_ptr<cartesian> currentForce {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldForce {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> currentVelocity {std::make_unique<cartesian>(cartesian::Zero())};
    std::unique_ptr<cartesian> oldVelocity {std::make_unique<cartesian>(cartesian::Zero())};

    // lock non const member access.
    // TO BE QUESTIONED
    tbb::spin_mutex mutex {};

private:
};