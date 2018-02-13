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

#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <tbb/spin_rw_mutex.h>
#include <tbb/scalable_allocator.h>


namespace enhance
{
    template
    <
        typename T, 
        template <typename, typename = tbb::scalable_allocator<T>> class Container
    >
    class ConcurrentContainer
    {
    public:
        typedef T MemberType;
        typedef Container<T> MemberList;
        typedef tbb::spin_rw_mutex mutex_type;
        typedef typename MemberList::iterator iterator;
        typedef typename MemberList::const_iterator const_iterator;

        ConcurrentContainer<T,Container>& operator=(ConcurrentContainer<T,Container>&&);

        // member management
        template<typename W>
        void add(W&&);

        template<typename W>
        void add_if_new(W&&);

        template<typename W>
        bool contains(const W&);
        
        template<typename W>
        void remove(W&);

        // information
        std::size_t size() const;

        // iteration
        typename MemberList::iterator begin();
        typename MemberList::iterator end();
        typename MemberList::const_iterator begin() const;
        typename MemberList::const_iterator end() const;

    protected:
        MemberList members;

    private:
        // thread safety
        mutex_type mutex;
    };



    template<typename T>
    using ConcurrentVector = ConcurrentContainer<T,std::vector>;
    template<typename T>
    using ConcurrentDeque = ConcurrentContainer<T,std::deque>;
    template<typename T>
    using ConcurrentList = ConcurrentContainer<T,std::list>;
    


    template<typename T,template <typename, typename> class Container>
    inline ConcurrentContainer<T,Container>& ConcurrentContainer<T,Container>::operator=(ConcurrentContainer<T,Container>&& other)
    {
        if(this!=&other)
        {
            members = other.members;
            other.members.clear();
        }
        return *this;
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    void ConcurrentContainer<T,Container>::add(W&& other)
    {
        mutex_type::scoped_lock lock(mutex, true);
        members.emplace_back(other);
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    void ConcurrentContainer<T,Container>::add_if_new(W&& other)
    {
        mutex_type::scoped_lock lock(mutex, false);
        if( !std::any_of(begin(),end(),[&](const T& member){ return member == other; }) ) 
        {
            lock.upgrade_to_writer();
            members.emplace_back(other);
        }
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    bool ConcurrentContainer<T,Container>::contains(const W& other)
    {
        mutex_type::scoped_lock lock(mutex, false);
        return std::any_of(begin(),end(),[&](const T& member){ return member == other; });
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    void ConcurrentContainer<T,Container>::remove(W& other)
    {
        mutex_type::scoped_lock lock(mutex, true);
        members.erase(std::remove(begin(), end(), other), end());
    }



    template<typename T,template <typename, typename> class Container>
    inline std::size_t ConcurrentContainer<T,Container>::size() const
    {
        return members.size();
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::MemberList::iterator ConcurrentContainer<T,Container>::begin()
    {
        return std::begin(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::MemberList::iterator ConcurrentContainer<T,Container>::end()
    {
        return std::end(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::MemberList::const_iterator ConcurrentContainer<T,Container>::begin() const
    {
        return std::cbegin(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::MemberList::const_iterator ConcurrentContainer<T,Container>::end() const
    {
        return std::cend(members);
    }
}