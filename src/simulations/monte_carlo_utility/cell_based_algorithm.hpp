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

#include "enhance/random.hpp"
#include "enhance/iterator_utility.hpp"
#include "enhance/parallel.hpp"
#include "cell_container.hpp"
#include <iterator>
#include <thread>
#include <chrono>
#include <tbb/parallel_for_each.h>
#include <tbb/scalable_allocator.h>
// #include <tbb/threads.h>



class CellBasedAlgorithm
{
public:
    template<class CELLCONTAINER>
    static void preparation( CELLCONTAINER& cells)
    {
        tbb::parallel_for_each(std::begin(cells), std::end(cells), [](typename CELLCONTAINER::cell_type& cell)
        {
            cell.state = CellState::IDLE;
        });
    }


    
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
        #ifndef NDEBUG
        using namespace std::chrono_literals;
        std::this_thread::sleep_for(100ms);
        #endif
        assert( cells.template allInState<CellState::STATE::FINISHED>() );
    }



    template<class CELLCONTAINER>
    static void reorder( CELLCONTAINER& cells)
    { 
        vesDEBUG(__PRETTY_FUNCTION__)
        typedef typename CELLCONTAINER::cell_type cell_type;
        typedef typename cell_type::particle_type particle_type;

        tbb::parallel_for_each(std::begin(cells), std::end(cells), [](cell_type& cell)
        {
            auto leavers = cell.getLeavers();
            
            for(particle_type& leaver : leavers)
            {
            // #ifndef NDEBUG
                bool was_added = false;
                for(cell_type& proximity_cell : cell.getProximity())
                {
                    was_added = proximity_cell.try_add(leaver);
                    if(was_added) break;
                }
                assert(was_added);
                // this somehow doesn't work. don't know why
            // #else
            //     for(cell_type& proximity_cell : cell.getProximity())
            //     {
            //         if(proximity_cell.contains(leaver.coords()))
            //         {
            //             proximity_cell.try_add(leaver);
            //         }
            //     }
            // #endif
                cell.removeParticle(leaver);
                assert(!cell.contains(&leaver));
            }
        });
    }
};