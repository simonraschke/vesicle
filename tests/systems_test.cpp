#define BOOST_TEST_MODULE systems
#include <boost/test/included/unit_test.hpp>
#include "systems/controller.hpp"


BOOST_AUTO_TEST_SUITE(systems)


BOOST_AUTO_TEST_CASE(periodic_box)
{
    Box<PERIODIC::ON> periodic;
    periodic.setParameters(Parameters());
    periodic.setLengthX(10);
    periodic.setLengthY(10);
    periodic.setLengthZ(10);

    Eigen::Vector3f a (1,1,1);
    Eigen::Vector3f b (9,1,1);

    BOOST_CHECK_MESSAGE(periodic.contains(a), "contains a");
    BOOST_CHECK_MESSAGE(periodic.contains(b), "contains b");

    // the next two should be true, because Box::contains calls Box::scaleDown
    BOOST_CHECK_MESSAGE(periodic.contains(-a), "not contains -a");
    BOOST_CHECK_MESSAGE(periodic.contains(-b), "not contains -b");
    
    BOOST_CHECK_MESSAGE(periodic.distance(a,b)-2 < 1e-3, "distance");
}



BOOST_AUTO_TEST_CASE(normal_box)
{
    Box<PERIODIC::OFF> normal;
    normal.setParameters(Parameters());
    normal.setLengthX(10);
    normal.setLengthY(10);
    normal.setLengthZ(10);

    Eigen::Vector3f a (1,1,1);
    Eigen::Vector3f b (9,1,1);

    BOOST_CHECK_MESSAGE(normal.contains(a), "contains a");
    BOOST_CHECK_MESSAGE(normal.contains(b), "contains b");
    BOOST_CHECK_MESSAGE(!normal.contains(-a), "not contains -a");
    BOOST_CHECK_MESSAGE(!normal.contains(-b), "not contains -b");
    BOOST_CHECK_MESSAGE(normal.distance(a,b)-8 < 1e-3, "distance");
}



BOOST_AUTO_TEST_CASE(system)
{
    System system;

    BOOST_CHECK_MESSAGE(system.getTime() < 1e-3, "time before");
    system.addTime(0.5);
    BOOST_CHECK_MESSAGE(system.getTime()-0.5 < 1e-3, "time after");
}



BOOST_AUTO_TEST_CASE(system_distribute_from_gro)
{
    const char* argv[3] = {nullptr,"--config","../../tests/test_config.ini"};

    System system;
    {
        Parameters prms;
        prms.programOptions.read(3,argv);
        prms.setup();
        system.setParameters(prms);
    }
    
    system.addParticles(ParticleFactory<ParticleMobile>(system.getParameters().mobile));
    TrajectoryDistributor dist;
    dist.setParameters(system.getParameters());
    // dist.setPath("../../tests/test_config.ini");
    dist(&system.getParticles());

    BOOST_CHECK_MESSAGE( system.getParticles().size() == 2, "number of particles "+std::to_string(system.getParticles().size()) ); 
}


BOOST_AUTO_TEST_SUITE_END()