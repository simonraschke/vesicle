#define likely(x)   __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)


#include <vector>
#include <memory>
#include "particles/particle.hpp"
#define PARTICLERANGE std::vector<std::unique_ptr<Particle>>


#include <eigen3/Eigen/Core>
//-------------------Eigen::IOFormat( prec, flag,                 coeffSep, rowSep, rowPre, rowSuf, matPre, matSuf )
#define ROWFORMAT    Eigen::IOFormat( 4,    Eigen::DontAlignCols, ", ",     ", ",   "",      "",    " ",    " " )
#define PYTHONFORMAT Eigen::IOFormat( 4,    0,                    ", ",     "\n",   "[",     "]",   "[",    "]" )