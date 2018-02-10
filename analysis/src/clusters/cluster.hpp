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

#include "enhance/observer_ptr.hpp"
#include "particles/particle_simple.hpp"
#include <deque>
#include <utility> // pair
#include <algorithm>
#include <tbb/spin_mutex.h>



// container holding pointers to particles
// used for DBSCAN
// thread safe interface
class Cluster
{
public:
    // typedef enhance::observer_ptr<Particle> Member_type;
    // typedef std::vector<Member_type> MemberList;
    // typedef std::pair<Member_type,MemberList> MemberNeighbourPair;
    // typedef std::vector<MemberNeighbourPair> MembersNeighbourList;
    typedef ParticleSimple MemberType;
    typedef std::deque<ParticleSimple> MemberList;
    typedef tbb::spin_mutex mutex_type;

    // member management
    void add(Particle&);
    bool contains(const Particle*) const;
    bool contains(const ParticleSimple&) const;
    void add_if_new(Particle&);
    // void expandNeighbourhood(MemberNeighbourPair&, Particle&);
    void remove(Particle*);
    void remove(const MemberType&);

    // preparation for cluster algorithms
    template<typename FUNCTOR>
    void constructNeighbourhood(PARTICLERANGE&, FUNCTOR&&);

    // iteration
    MemberList::iterator begin();
    MemberList::iterator end();
    MemberList::const_iterator begin() const;
    MemberList::const_iterator end() const;

protected:
private:
    // data
    // MembersNeighbourList member_neighbour_map;
    MemberList members;

    // thread safety
    mutex_type mutex;
};



// template<typename FUNCTOR>
// void Cluster::constructNeighbourhood(PARTICLERANGE& range, FUNCTOR&& func)
// {
    // for(auto& pair : member_neighbour_map)
    // {
    //     for(const auto& target : range)
    //     {
    //         if(func(pair.first,target))
    //         {
    //             expandNeighbourhood(pair, target.get());
    //         }
    //     }
    // }
// }