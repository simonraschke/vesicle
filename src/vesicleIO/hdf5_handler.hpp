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
#include "enhance/output_utility.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <memory>
#include <H5Cpp.h>

#include <variant>
#include <typeindex>
#include <type_traits>





class HDF5Handler
    : public ParameterDependentComponent
{
public:

    ~HDF5Handler();

    void setFileName(H5std_string);

    template<typename...Args>
    void createDataset(H5::PredType, H5std_string, Args...);

    void writeToDataset(H5std_string);

    bool isSupported(H5::DataSet);
    bool isSupported(H5::DataType);
    

protected:
    std::unique_ptr<H5std_string> file_name {nullptr};
    std::unique_ptr<H5::H5File> FILE {nullptr};
};





template<typename...Args>
inline void HDF5Handler::createDataset(H5::PredType T, H5std_string dataset_name, Args...ndims)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    static_assert(enhance::variadic::size<Args...>::value == 2, "enhance::SizeOfParameterPack(ndims...) != DIMENSIONS");
    static_assert(enhance::variadic::all_type<int,Args...>::value, "parameter pack does not hold all type int");
    assert(file_name);
    assert(FILE);

    // if(!supported_types.count(std::type_index(typeid(T))))
    //     throw std::logic_error(enhance::toStringViaStream("type" , typeid(T).name(), " not supported"));
    // else
    // {
    //     vesDEBUG(enhance::toStringViaStream("type" , typeid(T).name(), " supported"))
    // }

    hsize_t dimension_sizes[2] = {static_cast<hsize_t>(ndims)...};
	H5::DataSpace dataspace(2, dimension_sizes);
    H5::DataSet dataset = FILE->createDataSet(dataset_name, T, dataspace);

    if(!isSupported(dataset))
        throw std::logic_error("type of dataset not supported");
}