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

#include "systems/box.hpp"
#include "enhance/observer_ptr.hpp"
// #include <array>
#include <vector>
#include <memory>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>



template<typename T>
class Cell
    // : public Box<PERIODIC::ON>
{
public:
    typedef T particle_type;

protected:
    std::vector<enhance::observer_ptr<particle_type>> members{};
    std::vector<enhance::observer_ptr<const Cell<particle_type>>> proximity {};
    std::vector<enhance::observer_ptr<const Cell<particle_type>>> region {};

public:
    // bounding box access
    void setBoundaries(const Eigen::Vector3f&, const Eigen::Vector3f&);
    const Eigen::AlignedBox<float,3>& getBoundaries() const;

    bool contains(const Eigen::Vector3f&) const;
    bool contains(const std::unique_ptr<particle_type>&) const;
    bool contains(const enhance::observer_ptr<particle_type>&) const;

    // template<typename P>
    // bool isNeighbour(const Cell<P>&) const;

    // proximity and region access
    template<typename CRITERION>
    void setupProximityAndRegion(const std::vector<Cell<particle_type>>&, CRITERION&&);
    inline const decltype(proximity)& getProximity() const { return proximity; }
    inline const decltype(region)& getRegion() const { return region; }

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
    return std::any_of(std::cbegin(members),std::cend(members), p);
}



template<typename T>
inline bool Cell<T>::contains(const enhance::observer_ptr<particle_type>& p) const
{
    return std::any_of(std::cbegin(members),std::cend(members), p);
}



template<typename T>
template<typename CRITERION>
inline void Cell<T>::setupProximityAndRegion(const std::vector<Cell>& cells, CRITERION&& criterion)
{
    for(const auto& cell : cells)
    {
        if(std::addressof(cell) == std::addressof(*this))
            region.emplace_back( enhance::make_observer<const Cell<particle_type>>(&cell) );

        else if(criterion(*this,cell))
        {
            proximity.emplace_back( enhance::make_observer<const Cell<particle_type>>(&cell) );
            region.emplace_back( enhance::make_observer<const Cell<particle_type>>(&cell) );
        }
    }
}



// template<typename T>
// template<typename P>
// inline bool Cell<T>::isNeighbour(const Cell<P>& other) const
// {
//     Eigen::Vector3f connection_vector(bounding_box.center() - other.getBoundaries().center());
//     return true;
// }