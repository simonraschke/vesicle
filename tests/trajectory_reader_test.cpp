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

    {
        auto frame_lines = reader.getFrame(-1).second;
        auto it = std::begin(frame_lines);

        BOOST_CHECK_MESSAGE(frame_lines.size() == 7, "reader should have 7 lines, but has " + std::to_string(frame_lines.size()));

        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "FRAMEBEGIN t=0.0170"), *it );

        std::advance(it,1);
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "4"), *it );
        BOOST_CHECK_MESSAGE(!boost::algorithm::contains(*it, "."), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "1MOBIL    A    1   1.500   1.500   1.500  0.1487  0.0505 -0.1083"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "1MOBIL    B    2   0.500   0.500   0.500 -0.1487 -0.0505  0.1083"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "2MOBIL    A    3  -1.500  -1.500  -1.500  0.0721  0.0756  0.1707"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "2MOBIL    B    4  -0.500  -0.500  -0.500 -0.0721 -0.0756 -0.1707"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "10 10 10"), *it );
    }
    {
        // Any character besides zero must occur at least once
        BOOST_CHECK(reader.getMatches(std::regex("[^0]+")).count(1));
        auto frame_lines = reader.getMatches(std::regex("[^0]+")).rbegin()->second;
        auto it = std::begin(frame_lines);

        BOOST_CHECK_MESSAGE(frame_lines.size() == 7, "reader should have 7 lines, but has " + std::to_string(frame_lines.size()));

        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "FRAMEBEGIN t=0.0170"), *it );

        std::advance(it,1);
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "4"), *it );
        BOOST_CHECK_MESSAGE(!boost::algorithm::contains(*it, "."), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "1MOBIL    A    1   1.500   1.500   1.500  0.1487  0.0505 -0.1083"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "1MOBIL    B    2   0.500   0.500   0.500 -0.1487 -0.0505  0.1083"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "2MOBIL    A    3  -1.500  -1.500  -1.500  0.0721  0.0756  0.1707"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "2MOBIL    B    4  -0.500  -0.500  -0.500 -0.0721 -0.0756 -0.1707"), *it );
        
        std::advance(it,1);
        BOOST_CHECK(it != std::end(frame_lines));
        BOOST_CHECK_MESSAGE(boost::algorithm::contains(*it, "10 10 10"), *it );
    }
}


BOOST_AUTO_TEST_SUITE_END()