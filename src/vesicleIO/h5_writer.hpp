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

#include <boost/multi_array.hpp>

#ifdef __clang_major__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexceptions"
#pragma clang diagnostic ignored "-Wunused-parameter"
#include "h5xx/h5xx.hpp"
#pragma clang diagnostic pop
#elif  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wexceptions"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "h5xx/h5xx.hpp"
#pragma GCC diagnostic pop
#else
    #error no valid compiler
#endif

#include "trajectory.hpp"



struct TrajectoryWriterH5
    : public TrajectoryWriter
{   
    typedef float real;
    typedef boost::multi_array<std::uint16_t,1> array1d_t;
    typedef boost::multi_array<real,2> array2d_t;

    TrajectoryWriterH5();

    virtual void setPath(PATH) override;
    virtual void write(const float&, bool=false) override;  
    virtual void setAnisotropic(bool) override; 

protected:
    virtual void makeStartFileVMD() const override;

    array1d_t getResidueTypes() const;
    array2d_t getPositions() const;
    array2d_t getOrientations() const;

    h5xx::file h5file;
};