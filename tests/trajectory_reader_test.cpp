#define BOOST_TEST_MODULE trajectory_reader
#include <boost/test/included/unit_test.hpp>
#include "vesicleIO/gro_reader.hpp"


BOOST_AUTO_TEST_SUITE(trajectory_reader)


BOOST_AUTO_TEST_CASE(reader_gro)
{
    TrajectoryReaderGro reader;
    reader.setPath("../../tests/test_trajectory.gro");
    BOOST_CHECK_MESSAGE( reader.isOpen(), reader.getFilePath());
    reader.readAllFrames();

    auto frame_lines = reader.getFrame().second;
    auto it = std::begin(frame_lines);

    BOOST_CHECK_MESSAGE(frame_lines.size() == 7, "reader should have 7 lines, but has " + std::to_string(frame_lines.size()));

    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "FRAMEBEGIN t=0.0170"), *it );

    std::advance(it,1);
    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "4"), *it );
    BOOST_CHECK_MESSAGE(!boost::algorithm::contains(*it, "."), *it );
    
    std::advance(it,1);
    BOOST_CHECK(it != std::end(frame_lines));
    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "1MOBIL    A    1   0.655  -0.445  -0.761  0.1487  0.0505 -0.1083"), *it );
    
    std::advance(it,1);
    BOOST_CHECK(it != std::end(frame_lines));
    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "1MOBIL    B    2   1.352  -0.797  -0.136 -0.1487 -0.0505  0.1083"), *it );
    
    std::advance(it,1);
    BOOST_CHECK(it != std::end(frame_lines));
    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "2MOBIL    A    3  -0.679  -0.067   1.380  0.0721  0.0756  0.1707"), *it );
    
    std::advance(it,1);
    BOOST_CHECK(it != std::end(frame_lines));
    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "2MOBIL    B    4  -0.670   0.885   1.072 -0.0721 -0.0756 -0.1707"), *it );
    
    std::advance(it,1);
    BOOST_CHECK(it != std::end(frame_lines));
    BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "3.4200 3.4200 3.4200"), *it );
}


BOOST_AUTO_TEST_SUITE_END()