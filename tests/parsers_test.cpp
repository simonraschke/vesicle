#define BOOST_TEST_MODULE parsers
#include <boost/test/included/unit_test.hpp>
#include "systems/system.hpp"
// #include "parsers/potential_energy.hpp"
#include "parsers/surface_reconstruction.hpp"


BOOST_AUTO_TEST_SUITE(parsers)



BOOST_AUTO_TEST_CASE(surface_reconstruction_test)
{
    const char* argv[3] = {nullptr,"--config","../../tests/test_config.ini"};

    System system;
    Parameters prms;
    prms.read(3,argv);
    prms.setup();
    prms.mobile = 10;
    prms.analysis_cluster_volume_extension = -1;
    system.setParameters(prms);
    
    system.addParticles(ParticleFactory<ParticleMobile>(system.getParameters().mobile));
    system.getParticles()[0]->setCoords(Eigen::Vector3f(1,1,1));
    system.getParticles()[1]->setCoords(Eigen::Vector3f(2,1,1));
    system.getParticles()[2]->setCoords(Eigen::Vector3f(1,2,1));
    system.getParticles()[3]->setCoords(Eigen::Vector3f(1,1,2));
    system.getParticles()[4]->setCoords(Eigen::Vector3f(2,2,1));
    system.getParticles()[5]->setCoords(Eigen::Vector3f(1,2,2));
    system.getParticles()[6]->setCoords(Eigen::Vector3f(2,1,2));
    system.getParticles()[7]->setCoords(Eigen::Vector3f(2,2,2));
    system.getParticles()[8]->setCoords(Eigen::Vector3f(2,1,3));
    system.getParticles()[9]->setCoords(Eigen::Vector3f(2,2,3));
    
    // ClusterParser<PERIODIC::ON> clusters;
    // clusters.setParameters(system.getParameters());
    // clusters.setTarget(system.getParticles().begin(), system.getParticles().end());
    // clusters.DBSCANrecursive(1,3);

    // ClusterStructureParser structure(*clusters.begin());
    // structure.setParameters(system.getParameters());
    // structure.parse();
    // structure.printXML("../../tests/largest.vtu");
    
    // auto centre = structure.getCenter();
    // BOOST_CHECK_MESSAGE( std::abs(centre(0) - 1.6) < 1e-3, "centre x calculation wrong. " + std::to_string(centre(0)) );
    // BOOST_CHECK_MESSAGE( std::abs(centre(1) - 1.5) < 1e-3, "centre y calculation wrong. " + std::to_string(centre(1)) );
    // BOOST_CHECK_MESSAGE( std::abs(centre(2) - 1.8) < 1e-3, "centre z calculation wrong. " + std::to_string(centre(2)) );

    // BOOST_CHECK_MESSAGE( structure.getNumMembers() == 10, "members of structure wrong: " + std::to_string(structure.getNumMembers()) );
    // BOOST_CHECK_MESSAGE( std::abs(structure.getVolume() - 1.5) < 1e-3, "volume of structure wrong: " + std::to_string(structure.getVolume()) );
    // BOOST_CHECK_MESSAGE( std::abs(structure.getSurfaceArea() - 8.414) < 1e-3, "surface area of structure wrong: " + std::to_string(structure.getSurfaceArea()) );
}



BOOST_AUTO_TEST_SUITE_END()