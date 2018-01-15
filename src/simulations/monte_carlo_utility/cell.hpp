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
#include "cell_state.hpp"
#include "systems/box.hpp"
#include <vector>
#include <memory>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <tbb/spin_rw_mutex.h>



template<typename T>
class Cell
    : public Box<PERIODIC::ON>
{
public:
    typedef T particle_type;

protected:
    std::vector<std::reference_wrapper<particle_type>> particles{};
    std::vector<std::reference_wrapper<Cell<particle_type>>> proximity {};
    std::vector<std::reference_wrapper<Cell<particle_type>>> region {};

public:
    // bounding box access
    void setBoundaries(const Eigen::Vector3f&, const Eigen::Vector3f&);
    const Eigen::AlignedBox<float,3>& getBoundaries() const;

    bool contains(const typename particle_type::cartesian&) const;
    bool contains(const particle_type*) ;

    // proximity and region access
    template<typename CONTAINER, typename CRITERION>
    void setupProximityAndRegion(CONTAINER&, CRITERION&&);
    inline const decltype(proximity)& getProximity() const { return proximity; }
    inline const decltype(region)& getRegion() const { return region; }

    // member access
    void clearParticles();
    void removeParticle(const particle_type&);
    bool try_add(particle_type&);
    std::vector<std::reference_wrapper<particle_type>> getLeavers();
    std::size_t size() const { return particles.size(); }
    constexpr auto begin() const {return std::begin(particles); }
    constexpr auto end() const {return std::end(particles); }
    constexpr auto cbegin() const {return std::cbegin(particles); }
    constexpr auto cend() const {return std::cend(particles); }

    //state
    CellState state {};

    template<CellState::STATE S>
    bool proximityAllInState() const;
    template<CellState::STATE S>
    bool proximityNoneInState() const;

    template<CellState::STATE S>
    bool regionAllInState() const;
    template<CellState::STATE S>
    bool regionNoneInState() const;

private:
    Eigen::AlignedBox<float,3> bounding_box {};
    tbb::spin_rw_mutex particles_access_mutex {};
};



template<typename T>
inline void Cell<T>::setBoundaries(const Eigen::Vector3f& from, const Eigen::Vector3f& to)
{
    bounding_box.setEmpty();
    bounding_box.extend(from);
    bounding_box.extend(to);
}



template<typename T>
inline bool Cell<T>::contains(const typename particle_type::cartesian& coords) const
{
    // vesDEBUG( " scaled coords " << scaleDown(coords).format(ROWFORMAT) )
    return bounding_box.contains(scaleDown(coords));
}



template<typename T>
inline const Eigen::AlignedBox<float,3>& Cell<T>::getBoundaries() const
{
    return bounding_box;
}



template<typename T>
inline bool Cell<T>::contains(const particle_type* other) 
{
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, false);
    assert(other);
    return std::any_of(std::cbegin(particles),std::cend(particles), [&](const particle_type& p ){ return other == &p; });
}



template<typename T>
template<typename CONTAINER, typename CRITERION>
inline void Cell<T>::setupProximityAndRegion(CONTAINER& cells, CRITERION&& criterion)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    for(Cell& cell : cells)
    {
        if(std::addressof(cell) == std::addressof(*this))
            region.emplace_back( std::ref(cell) );
        else if(criterion(*this,cell))
        {
            proximity.emplace_back( std::ref(cell) );
            region.emplace_back( std::ref(cell) );
        }
    }
}



template<typename T>
inline void Cell<T>::clearParticles()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
    particles.clear();
}



template<typename T>
inline void Cell<T>::removeParticle(const particle_type& to_remove)
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
    particles.erase( std::remove_if(std::begin(particles), std::end(particles), [&](const particle_type& to_compare)
    { 
        return std::addressof(to_remove) == std::addressof(to_compare);
    ;}), particles.end() );
}



template<typename T>
inline bool Cell<T>::try_add(particle_type& particle)
{
    if(!contains(&particle) && contains(particle.coords()))
    {   
        tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
        particles.emplace_back(std::ref(particle));
        lock.release();
        // vesDEBUG(__func__ << "  Cell contains particle after insertion " << std::boolalpha << contains(&particle))
        assert(contains(&particle));
        return true;
    }
    else
    {
        // vesDEBUG(__func__ << "  Cell contains particle after NO insertion " << std::boolalpha << contains(&particle))
        assert(!contains(&particle));
        return false;
    }
}



template<typename T>
inline auto Cell<T>::getLeavers() -> std::vector<std::reference_wrapper<particle_type>> 
{
    tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, false);
    decltype(particles) leavers;
    for(particle_type& particle : *this)
    {
        if(!contains(particle.coords()))
        {
            leavers.emplace_back(std::ref(particle));
            assert(!contains(leavers.back().get().coords()));
        }
    }
    return leavers;
}



template<typename T>
template<CellState::STATE S>
bool Cell<T>::proximityAllInState() const
{
    return std::all_of(std::begin(proximity),std::end(proximity), [](const Cell<T>& cell){ return cell.state == S; } );
}



template<typename T>
template<CellState::STATE S>
bool Cell<T>::proximityNoneInState() const
{
    return std::none_of(std::begin(proximity),std::end(proximity), [](const Cell<T>& cell){ return cell.state == S; } );
}



template<typename T>
template<CellState::STATE S>
bool Cell<T>::regionAllInState() const
{
    return std::all_of(std::begin(region),std::end(region), [](const Cell<T>& cell){ return cell.state == S; } );
}



template<typename T>
template<CellState::STATE S>
bool Cell<T>::regionNoneInState() const
{
    return std::none_of(std::begin(region),std::end(region), [](const Cell<T>& cell){ return cell.state == S; } );
}