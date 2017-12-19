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

#include "enhance/random.hpp"
#include "enhance/iterator_utility.hpp"
#include "enhance/parallel.hpp"
#include "cell_container.hpp"
#include <iterator>



class CellBasedAlgorithm
{
public:
    template<class CELLCONTAINER, typename FUNCTOR>
    static void step( CELLCONTAINER& cells, FUNCTOR&& func )
    {
        vesDEBUG(__PRETTY_FUNCTION__)
        static enhance::scoped_root_dummy ROOT;
        while( ! (cells.template allInState<CellState::FINISHED>()) )
        {
            enhance::for_each_from(
                std::begin(cells), 
                std::end(cells), 
                std::begin(cells) + enhance::random<std::size_t>(0,std::distance(std::begin(cells),std::end(cells))-1), 
                [&](typename CELLCONTAINER::cell_type& cell)
            {
                if( 
                    cell.template regionNoneInState<CellState::BLOCKED>() && 
                    cell.state == CellState::IDLE
                )
                {
                    cell.state = CellState::BLOCKED;
                    
                    assert( cell.state == CellState::BLOCKED );
                    assert( cell.template proximityNoneInState<CellState::BLOCKED>() );
                    
                    ROOT.enqueue_child( [&]
                    {
                        assert( cell.state == CellState::BLOCKED );
                        func( cell ); 
                        cell.state = CellState::FINISHED;
                        assert( cell.state == CellState::FINISHED );
                    } );
                }
            });
        }
        assert( cells.template allInState<CellState::STATE::FINISHED>() );
    }
};