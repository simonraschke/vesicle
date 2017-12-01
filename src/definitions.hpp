#pragma once

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

#include <csignal>

#ifndef NDEBUG
    #include <iostream>
    #define vesDEBUG(x) {std::cerr << "[DEBUG] "; do { std::cerr << x; } while (0); std::cerr << '\n';}
#else
    #define vesDEBUG(x)
#endif
#define vesLOG(x) {std::clog << "[LOG] "; do { std::clog << x; } while (0); std::clog << '\n';}
#define vesWARNING(x) {std::clog << "[WARNING] "; do { std::clog << x; } while (0); std::clog << '\n';}
#define vesCRITICAL(x) {std::cerr << "[ERROR] "<< __FILE__ <<":" << __LINE__ << "  "; do { std::cerr << x; } while (0); std::cerr <<" raising SIGABRT\n"; std::raise(SIGABRT);}


#include <vector>
#include <memory>
#include "particles/particle.hpp"
#define PARTICLERANGE std::vector<std::unique_ptr<Particle>>


#include <eigen3/Eigen/Core>
//-------------------Eigen::IOFormat( prec, flag,                 coeffSep, rowSep, rowPre, rowSuf, matPre, matSuf )
#define ROWFORMAT    Eigen::IOFormat( 3,    Eigen::DontAlignCols, ", ",     " ",    " ",    "",     " ",    " " )
#define PYTHONFORMAT Eigen::IOFormat( 4,    0,                    ", ",     "\n",   "[",    "]",    "[",    "]" )


#include <tbb/parallel_reduce.h>
#define PARALLEL_REDUCE(type, cont, functor)  tbb::parallel_reduce \
    (tbb::blocked_range<decltype(cont)::const_iterator>( std::cbegin(cont), std::cend(cont) ), \
    (type)0 , [&](auto& r, type i) \
    { \
        return i + std::accumulate(std::begin(r), std::end(r), (type)0, functor); \
    },\
    std::plus<type>())

