#include "hdf5_handler.hpp"




HDF5Handler::~HDF5Handler()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    FILE.reset(nullptr);
}



void HDF5Handler::setFileName(H5std_string name)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(!name.empty());
    file_name = std::make_unique<H5std_string>(name);
    FILE = std::make_unique<H5::H5File>(*file_name, H5F_ACC_TRUNC);
    assert(FILE);
}



bool HDF5Handler::isSupported(H5::DataSet dataset)
{
    return isSupported(dataset.getDataType());
}



bool HDF5Handler::isSupported(H5::DataType datatype)
{
    if(datatype == H5::PredType::NATIVE_FLOAT)
        return true;
    else if(datatype == H5::PredType::NATIVE_HSIZE)
        return true;
    else if(datatype == H5::PredType::NATIVE_INT32)
        return true;
    else if(datatype == H5::PredType::NATIVE_INT64)
        return true;
    else 
        return false;
}



void HDF5Handler::writeToDataset(H5std_string name)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    assert(!name.empty());

	H5::DataSet dataset = FILE->openDataSet(name);
    H5::DataType datatype = dataset.getDataType();

    std::variant
    <
        std::array<std::array<int,10000>,9999>,
        std::array<std::array<long,10000>,9999>,
        std::array<std::array<float,10000>,9999>,
        std::array<std::array<hsize_t,10000>,9999>
    >
    data;    

    if(datatype == H5::PredType::NATIVE_FLOAT)
        data = std::array<std::array<float,10000>,9999>();
    else if(datatype == H5::PredType::NATIVE_HSIZE)
        data = std::array<std::array<hsize_t,10000>,9999>();
    else if(datatype == H5::PredType::NATIVE_INT32)
        data = std::array<std::array<int,10000>,9999>();
    else if(datatype == H5::PredType::NATIVE_INT64)
        data = std::array<std::array<long,10000>,9999>();
    else
        throw std::logic_error("type of dataset not supported");

    std::visit([&](auto&& arg)
    {
        for (int j = 0; j < 9999; j++)
        {
            for (int i = 0; i < 10000; i++)
            {
                using T = std::decay_t<decltype(arg)>;
                if constexpr (std::is_same_v<T, std::array<std::array<int,10000>,9999>>)
                    arg[j][i] = 1.337;
                else if constexpr (std::is_same_v<T, std::array<std::array<long,10000>,9999>>)
                    arg[j][i] = 1337;
                else if constexpr (std::is_same_v<T, std::array<std::array<hsize_t,10000>,9999>>)
                    arg[j][i] = 1337;
                else if constexpr (std::is_same_v<T, std::array<std::array<float,10000>,9999>>)
                    arg[j][i] = 1337;
                else
                    static_assert(enhance::always_false<T>::value, "non-exhaustive visitor!");
            }
        }

        try
        {
            dataset.write(arg.data()->data(), dataset.getDataType());
        }

        // catch failure caused by the H5File operations
        catch(H5::FileIException error)
        {
            vesDEBUG("H5::FileIException")
            error.printError();
        }

        // catch failure caused by the DataSet operations
        catch(H5::DataSetIException error)
        {
            vesDEBUG("H5::DataSetIException")
            error.printError();
        }
    }, data);
}
