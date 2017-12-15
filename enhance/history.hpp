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
