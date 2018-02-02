#include "hdf5_handler.hpp"



void HDF5Handler::setFileName(H5std_string name)
{
    assert(!name.empty());
    file_name = std::make_unique<H5std_string>(name);
    FILE = std::make_unique<H5::H5File>(*file_name, H5F_ACC_TRUNC);
    assert(FILE);
}



void HDF5Handler::writeToDataset(H5std_string name)
{
    assert(!name.empty());

	H5::DataSet dataset = FILE->openDataSet(name);

    int i, j;
    float data[4][6];	    // buffer for data to write
    for (j = 0; j < 4; j++)
    for (i = 0; i < 6; i++)
        data[j][i] = i * 6 + j + 1;
        // data[j][i] = 1;

	dataset.write(data, dataset.getDataType());
}
