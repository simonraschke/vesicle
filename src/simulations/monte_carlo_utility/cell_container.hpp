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

#include "cell.hpp"
#include "systems/box.hpp"
#include "enhance/random.hpp"
#include <tbb/parallel_for_each.h>
#include <deque>



template<typename CELL_MEM_T>
class CellContainer
    : public Box<PERIODIC::OFF>
    , public virtual ParameterDependentComponent
{
public:
    typedef Cell<CELL_MEM_T> cell_type;

    void setup();
    void reorder();

    template<typename P1, typename P2>
    bool areNeighbourCells(const Cell<P1>&, const Cell<P2>&) const;

    template<typename CONTAINER>
    void deployParticles(const CONTAINER&);

    constexpr auto begin() {return std::begin(cells); }
    constexpr auto end() {return std::end(cells); }

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
    const Eigen::Vector3f connection_vector = periodic.distance_vector(first.getBoundaries().center(), second.getBoundaries().center()).cwiseAbs();
    
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

    const std::size_t cells_x = std::trunc( getLengthX() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( getLengthX() / (min_edge) )  ;
    const float x_edge = getLengthX()/cells_x;

    const std::size_t cells_y = std::trunc( getLengthY() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( getLengthY() / (min_edge) )  ;
    const float y_edge = getLengthY()/cells_y;

    const std::size_t cells_z = std::trunc( getLengthZ() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( getLengthZ() / (min_edge) )  ;
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
    vesDEBUG(__PRETTY_FUNCTION__)
    return std::all_of(std::begin(cells),std::end(cells), [](const Cell<CELL_MEM_T>& cell){ return cell.state == S; } );
}



template<typename CELL_MEM_T>
template<CellState::STATE S>
bool CellContainer<CELL_MEM_T>::noneInState() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    return std::none_of(std::begin(cells),std::end(cells), [](const Cell<CELL_MEM_T>& cell){ return cell.state == S; } );
}



template<typename CELL_MEM_T>
void CellContainer<CELL_MEM_T>::reorder() 
{
    // this one MIGHT be not thread-safe at "put(L)"
    // should be now. spin_mutex in member interface

    tbb::parallel_for_each(std::begin(cells), std::end(cells), [&](Cell<CELL_MEM_T> & cell)
    {
        cell.clearParticles();
    });
//     tbb::parallel_for(0,(int)size_,1,[&](const SIZE& __ID)
//     {
//         module::Cell& CELL = cells_[__ID];
//         if(CELL.need_reorder==true)
//         {
//             scalable_vector<SIZE> leavers;
//             for(const auto& M : CELL.members())
//             {
//                 if(particles_[M]->left_cell) 
//                 {
//                     leavers.emplace_back(M);
//                 }
//             }
            
//             for(const auto& L : leavers)
//             {
//                 CELL.members.remove(L);
// //             }
// //             
// //             for(const auto& L : leavers)
// //             {
//                 put(L);
//                 particles_[L]->left_cell = false;
//             }
            
//             CELL.need_reorder = false;
//         }
//     });
}