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

#include "cell.hpp"
#include "systems/box.hpp"
#include "enhance/random.hpp"
#include <tbb/parallel_for_each.h>
#include <deque>


// class owning the cells and provide an interface for them
// is iterable
// template parameter will be the type of Cell
template<typename CELL_MEM_T>
class CellContainer
    : public Box<PERIODIC::OFF>
    , public virtual ParameterDependentComponent
{
public:
    typedef Cell<CELL_MEM_T> cell_type;

    // call before usage
    void setup();

    // check if two cells are neigbours
    template<typename P1, typename P2>
    bool areNeighbourCells(const Cell<P1>&, const Cell<P2>&) const;

    template<typename CONTAINER>
    void deployParticles(const CONTAINER&);

    constexpr auto begin() {return std::begin(cells); }
    constexpr auto end() {return std::end(cells); }

    std::size_t membersContained() const;

    // state
    template<CellState::STATE S>
    bool allInState() const;

    template<CellState::STATE S>
    bool noneInState() const;

private:
    std::deque<cell_type> cells {};
};



template<typename CELL_MEM_T>
template<typename P1, typename P2>
inline bool CellContainer<CELL_MEM_T>::areNeighbourCells(const Cell<P1>& first, const Cell<P2>& second) const
{
    Box<PERIODIC::ON> periodic;
    periodic.setParameters(getParameters());
    const Eigen::Vector3f connection_vector = periodic.distanceVector(first.getBoundaries().center(), second.getBoundaries().center()).cwiseAbs();
    
    if(std::addressof(first) == std::addressof(second))
    {
        vesDEBUG(__func__ << " " << connection_vector.format(ROWFORMAT) << " return:  false  (Reason: same Cell)")
        return false;
    }
    else if((connection_vector(0) < first.getBoundaries().sizes()(0) + 0.01) 
        &&  (connection_vector(1) < first.getBoundaries().sizes()(1) + 0.01) 
        &&  (connection_vector(2) < first.getBoundaries().sizes()(2) + 0.01) )
    { 
        vesDEBUG(__func__ << " " << connection_vector.format(ROWFORMAT) << " return:  true  (Reason: in reach)")
        return true;
    }
    else 
    { 
        vesDEBUG(__func__ << " " << connection_vector.format(ROWFORMAT) << " return:  false (Reason: not in reach)")
        return false;
    }
}



template<typename CELL_MEM_T>
inline void CellContainer<CELL_MEM_T>::setup()
{
    const float min_edge = getParameters().cell_min_edge;
    const std::size_t max_cells_dim = getParameters().max_cells_dim;

    vesDEBUG("box x y z: " << getLengthX()  << " " << getLengthY()  << " " << getLengthZ())
    vesDEBUG("minimum edge for a cell is " << min_edge)
    vesDEBUG("maximum number of cells per dimension is " << max_cells_dim)

    // calculate the amount of zells per dimension
    std::size_t x_helper = std::trunc( getLengthX() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( getLengthX() / (min_edge) )  ;
    std::size_t y_helper = std::trunc( getLengthY() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( getLengthY() / (min_edge) )  ;
    std::size_t z_helper = std::trunc( getLengthZ() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( getLengthZ() / (min_edge) )  ;

    // if dimension is smaller than min_edge it will be 0, therefor set to 1
    const std::size_t cells_x = x_helper == 0 ? 1 : x_helper;
    const std::size_t cells_y = y_helper == 0 ? 1 : y_helper;
    const std::size_t cells_z = z_helper == 0 ? 1 : z_helper;

    // edge of cell per dimension
    const float x_edge = getLengthX()/cells_x;
    const float y_edge = getLengthY()/cells_y;
    const float z_edge = getLengthZ()/cells_z;

    vesDEBUG("cells in x dimension: " << cells_x << " with edge: " << x_edge)
    vesDEBUG("cells in y dimension: " << cells_y << " with edge: " << y_edge)
    vesDEBUG("cells in z dimension: " << cells_z << " with edge: " << z_edge)
    vesDEBUG("cells overall: " << cells_x*cells_y*cells_z)

    for ( std::size_t x = 0; x < cells_x; ++x )
    for ( std::size_t y = 0; y < cells_y; ++y )
    for ( std::size_t z = 0; z < cells_z; ++z )
    {
        vesDEBUG("##### building cell: " << x << ' ' << y << ' ' << z)
        cells.emplace_back();
        cells.back().setParameters(getParameters());
        const Eigen::Vector3f from = Eigen::Vector3f( x_edge*x,     y_edge*y,     z_edge*z     ) - Eigen::Vector3f(0.01, 0.01, 0.01);
        const Eigen::Vector3f to   = Eigen::Vector3f( x_edge*(x+1), y_edge*(y+1), z_edge*(z+1) ) + Eigen::Vector3f(0.01, 0.01, 0.01);
        cells.back().setBoundaries(from,to);

        std::cerr.precision(3);
        std::cerr << std::fixed;
        vesDEBUG("  from:  " << from.format(ROWFORMAT))
        vesDEBUG("  to:    " << to.format(ROWFORMAT))
        vesDEBUG("  min    " << cells.back().getBoundaries().min().format(ROWFORMAT))
        vesDEBUG("  max    " << cells.back().getBoundaries().max().format(ROWFORMAT))
        vesDEBUG("  centre " << cells.back().getBoundaries().center().format(ROWFORMAT))
    }

    assert( cells.size() == cells_x*cells_y*cells_z );

    tbb::parallel_for_each(std::begin(cells), std::end(cells), [&](Cell<CELL_MEM_T> & cell)
    {
        cell.setupProximityAndRegion(cells, [&](const auto& a, const auto& b) { return areNeighbourCells(a,b); } );
        vesDEBUG(cell.getProximity().size());
        cell.state = CellState::IDLE;
    });

    #ifndef NDEBUG
    std::for_each(std::begin(cells), std::end(cells), [&](Cell<CELL_MEM_T> & cell)
    {
        vesDEBUG("size of cell proximity: " << cell.getProximity().size() )
        vesDEBUG("size of cell region:    " << cell.getRegion().size() )
        assert(cell.getProximity().size() == 26);
        assert(cell.getRegion().size() == 27);
        assert(cell.state == CellState::IDLE);
    });
    #endif
}



template<typename CELL_MEM_T>
std::size_t CellContainer<CELL_MEM_T>::membersContained() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return PARALLEL_REDUCE(std::size_t, cells, [](auto i, const cell_type& cell){ return i + cell.size(); });
}



template<typename CELL_MEM_T>
template<typename CONTAINER>
void CellContainer<CELL_MEM_T>::deployParticles(const CONTAINER& particles)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    tbb::parallel_for_each(std::begin(particles), std::end(particles), [&](const std::unique_ptr<CELL_MEM_T>& particle)
    {
        bool done = false;
        std::for_each(std::begin(cells), std::end(cells), [&](Cell<CELL_MEM_T> & cell)
        {
            if(done) return;
            done = cell.try_add(*particle);
        });
        assert(done);
    });
}



template<typename CELL_MEM_T>
template<CellState::STATE S>
bool CellContainer<CELL_MEM_T>::allInState() const
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    return std::all_of(std::begin(cells),std::end(cells), [](const Cell<CELL_MEM_T>& cell){ return cell.state == S; } );
}



template<typename CELL_MEM_T>
template<CellState::STATE S>
bool CellContainer<CELL_MEM_T>::noneInState() const
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    return std::none_of(std::begin(cells),std::end(cells), [](const Cell<CELL_MEM_T>& cell){ return cell.state == S; } );
}


