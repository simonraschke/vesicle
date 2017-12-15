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

// #include "particles/particle.hpp"
#include "enhance/observer_ptr.hpp"
#include <array>
#include <vector>
#include <memory>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>



template<typename T>
class Cell
{
public:
    typedef T particle_type;

    // bounding box access
    void setBoundaries(const Eigen::Vector3f&, const Eigen::Vector3f&);
    const Eigen::AlignedBox<float,3>& getBoundaries() const;

    bool contains(const Eigen::Vector3f&) const;
    bool contains(const std::unique_ptr<particle_type>&) const;
    bool contains(const enhance::observer_ptr<particle_type>&) const;

protected:
    std::vector<enhance::observer_ptr<particle_type>> members{};
    std::array<enhance::observer_ptr<Cell<particle_type>>,26> proximity {};

private:
    Eigen::AlignedBox<float,3> bounding_box {};
};



template<typename T>
inline void Cell<T>::setBoundaries(const Eigen::Vector3f& from, const Eigen::Vector3f& to)
{
    bounding_box.setEmpty();
    bounding_box.extend(from);
    bounding_box.extend(to);
}



template<typename T>
inline bool Cell<T>::contains(const Eigen::Vector3f& vec) const
{
    return bounding_box.contains(vec);
}



template<typename T>
inline const Eigen::AlignedBox<float,3>& Cell<T>::getBoundaries() const
{
    return bounding_box;
}



template<typename T>
inline bool Cell<T>::contains(const std::unique_ptr<particle_type>& p) const
{
    // return contains(p.coords());
    return std::any_of(std::cbegin(members),std::cend(members), p);
}



template<typename T>
inline bool Cell<T>::contains(const enhance::observer_ptr<particle_type>& p) const
{
    // return contains(p.coords());
    return std::any_of(std::cbegin(members),std::cend(members), p);
}