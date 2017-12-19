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

#include <memory>
#include <array>



namespace enhance
{

    template<std::size_t SIZE, typename T>
    class History
    {
    public:
        template<typename...Args>
        void append(Args&& ...);

    protected:
    private:
        std::size_t position {0};
        std::array<std::unique_ptr<T>, SIZE> storage {};
    };



    
    template<std::size_t SIZE, typename T>
    template<typename...Args>
    void History<SIZE,T>::append(Args&& ... args)
    {
        if(storage[position] == nullptr)
            storage[position] = std::make_unique<T>(std::forward<Args>(args)...);
        // else
            // *storage[position] = std::move(args);
    }

    // template<typename T>
    // template<typename...Args>
    // History::History(Args&& ... args)
    // {
        
    // }
}
