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

#include "parameters.hpp"
#include "enhance/compile_time_utility.hpp"
#include <iostream>
#include <string>
#include <memory>
#include <H5Cpp.h>




class HDF5Handler
    : public ParameterDependentComponent
{
public:
    void setFileName(H5std_string);

    template<hsize_t DIMENSIONS, typename...Args>
    void createDataset(H5::PredType, H5std_string, Args...);

    void writeToDataset(H5std_string);

protected:
    std::unique_ptr<H5std_string> file_name {nullptr};
    std::unique_ptr<H5::H5File> FILE {nullptr};
};





template<hsize_t DIMENSIONS, typename...Args>
inline void HDF5Handler::createDataset(H5::PredType T, H5std_string dataset_name, Args...ndims)
{
    static_assert(enhance::variadic::size<Args...>::value == DIMENSIONS, "enhance::SizeOfParameterPack(ndims...) != DIMENSIONS");
    static_assert(enhance::variadic::all_type<int,Args...>::value, "parameter pack does not hold all type int");
    assert(file_name);
    assert(FILE);

    hsize_t dimension_sizes[DIMENSIONS] = {static_cast<hsize_t>(ndims)...};
	H5::DataSpace dataspace( DIMENSIONS, dimension_sizes);
    H5::DataSet dataset = FILE->createDataSet(dataset_name, T, dataspace);
}