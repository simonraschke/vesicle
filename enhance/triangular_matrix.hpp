#pragma once

#include <vector>
#include <tbb/scalable_allocator.h>
#include <tbb/parallel_reduce.h>



namespace enhance
{
    template<typename T>
    class TriangularMatrix
    {
        std::size_t size_ = 0;
        std::vector<std::uint32_t,tbb::scalable_allocator<std::uint32_t>> sums_ {};
        std::vector<T,tbb::scalable_allocator<T>> data_ {};
        tbb::spin_mutex data_mutex;
        typedef tbb::spin_mutex::scoped_lock lock;
        
    public:
        explicit TriangularMatrix();
        TriangularMatrix(const std::size_t&);
        ~TriangularMatrix();

        T sum() noexcept;
        T minCoeff();
        T maxCoeff();
        void scale(const T&, const T&);
        void resize(const std::size_t&);

        T& operator[](const std::size_t&);
        const T& operator[](const std::size_t&) const;
        T& operator()(const std::size_t&, const std::size_t&);
        const T& operator()(const std::size_t&, const std::size_t&) const;
        
        auto begin() noexcept       { return data_.begin(); }
        auto cbegin() const noexcept{ return data_.cbegin(); }
        auto end() noexcept         { return data_.end(); }
        auto cend() const noexcept  { return data_.cend(); }
    };
}



template<typename T>
inline enhance::TriangularMatrix<T>::TriangularMatrix()
{    
}



template<typename T>
inline enhance::TriangularMatrix<T>::TriangularMatrix(const std::size_t& SIZE)
{    
    resize(SIZE);
}



template<typename T>
inline enhance::TriangularMatrix<T>::~TriangularMatrix()
{
}



template<typename T>
inline T& enhance::TriangularMatrix<T>::operator[](const std::size_t& i)
{
    return data_[i];
}



template<typename T>
inline const T& enhance::TriangularMatrix<T>::operator[](const std::size_t& i) const
{
    return data_[i];
}



template<typename T>
inline T& enhance::TriangularMatrix<T>::operator()(const std::size_t& _row, const std::size_t& _col)
{
    if(! (_col < _row)) return data_[sums_[_col]+_row];
    else return data_[sums_[_row]+_col];
}



template<typename T>
inline const T& enhance::TriangularMatrix<T>::operator()(const std::size_t& _row, const std::size_t& _col) const
{
    if(! (_col < _row)) return data_[sums_[_col]+_row];
    else return data_[sums_[_row]+_col];
}



template<typename T>
inline T enhance::TriangularMatrix<T>::minCoeff()
{
    lock(data_mutex);
    return *std::min_element(begin(), end());
}



template<typename T>
inline T enhance::TriangularMatrix<T>::maxCoeff()
{
    lock(data_mutex);
    return *std::max_element(begin(), end());
}



template<typename T>
inline void enhance::TriangularMatrix<T>::scale(const T& __from, const T& __to)
{
    lock(data_mutex);
    const T min = minCoeff();
    const T max = maxCoeff();
    const T scaling_factor = (__to - __from)/(max - min);
    
    for(std::size_t i = 0; i < size_; ++i)
    {
        data_[i] = data_[i]*scaling_factor + (__to + __from)/2;
    }
}



template<typename T>
inline void enhance::TriangularMatrix<T>::resize(const std::size_t& SIZE)
{
    lock(data_mutex);
    size_ = (SIZE*SIZE)/2;
    sums_.resize(SIZE);
    data_.resize(size_);

    sums_[0] = 0;
    for(std::uint32_t i(1); i < SIZE; ++i)
    {
        sums_[i] = sums_[i-1]+i;
    }
}



template<typename T>
inline T enhance::TriangularMatrix<T>::sum() noexcept
{
    lock(data_mutex);
    return tbb::parallel_reduce 
        (tbb::blocked_range<typename decltype(data_)::const_iterator>( std::cbegin(data_), std::cend(data_) ), 
        (T)0 , [&](auto& r, T i) 
        { 
            return i + std::accumulate(std::begin(r), std::end(r), (T)0, std::plus<T>()); 
        },
        std::plus<T>());
}
