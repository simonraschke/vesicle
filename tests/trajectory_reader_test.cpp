#define BOOST_TEST_MODULE trajectory_reader
#include <boost/test/included/unit_test.hpp>
#include "vesicleIO/gro_reader.hpp"


BOOST_AUTO_TEST_SUITE(trajectory_reader)


BOOST_AUTO_TEST_CASE(reader_gro)
{
    TrajectoryReaderGro reader;
    reader.setFilename("../tests/test_trajectory");
    reader.readLastFrame();
}


BOOST_AUTO_TEST_SUITE_END()