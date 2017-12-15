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



template<typename CELL_MEM_T>
class CellContainer
    : public Box<PERIODIC::OFF>
    , public virtual ParameterDependentComponent
{
public:
    void setup();

private:
    std::vector<Cell<CELL_MEM_T>> cells {};
};




template<typename CELL_MEM_T>
inline void CellContainer<CELL_MEM_T>::setup()
{
    float x_edge = 0.f;
    float y_edge = 0.f;
    float z_edge = 0.f;
    const float min_edge = getParameters().cell_min_edge;
    const std::size_t max_cells_dim = getParameters().max_cells_dim;

    if( min_edge/getLengthX() <= max_cells_dim)
        x_edge = min_edge;
    else
        x_edge = getLengthX()/max_cells_dim;
    std::size_t max_cells_x = std::round(getLengthX()/x_edge);

    if( min_edge/getLengthY() <= max_cells_dim)
        y_edge = min_edge;
    else
        y_edge = getLengthY()/max_cells_dim;
    std::size_t max_cells_y = std::round(getLengthY()/y_edge);

    if( min_edge/getLengthZ() <= max_cells_dim)
        z_edge = min_edge;
    else
        z_edge = getLengthZ()/max_cells_dim;
    std::size_t max_cells_z = std::round(getLengthZ()/z_edge);

    for ( std::size_t x = 0; x < max_cells_x; ++x )
    for ( std::size_t y = 0; y < max_cells_y; ++y )
    for ( std::size_t z = 0; z < max_cells_z; ++z )
    {
        cells.emplace_back( );
        auto from = Eigen::Vector3f( x_edge*x,     y_edge*y,     z_edge*z     ) - Eigen::Vector3f(0.01, 0.01, 0.01);
        auto to   = Eigen::Vector3f( x_edge*(x+1), y_edge*(y+1), z_edge*(z+1) ) + Eigen::Vector3f(0.01, 0.01, 0.01);
        cells.back().setBoundaries(from,to);
        #ifndef NDEBUG
            vesDEBUG("built cell:")
            vesDEBUG( "min    " << cells.back().getBoundaries().min().format(ROWFORMAT))
            vesDEBUG( "max    " << cells.back().getBoundaries().max().format(ROWFORMAT))
            vesDEBUG( "centre " << cells.back().getBoundaries().center().format(ROWFORMAT))
        #endif
    }

    std::terminate();
    // assert( this->size() == std::pow(cells_per_dim, 3) );
    
    // console::debug(size(), "cells were built");
    
    // setup_cell_environments();
    
    // #ifndef NDEBUG
    // parallel_for_each( [&](const cell_t& origin) { for_each([&](const cell_t& other) { if (origin == other) return; assert( origin != other ); } ); });
    // parallel_for_each( [&](const cell_t& origin) { origin.proximity.for_each( [&](const cell_t& other) { assert( origin != other ); }); });
    // parallel_for_each( [&](const cell_t& origin) { origin.region.for_each( [&](const cell_t& other)    { if (origin == other) return; assert( origin != other ); }); });
    // #endif
}