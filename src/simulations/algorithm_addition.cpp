#include "algorithm_addition.hpp"



void AlgorithmAddition::setup()
{
    vesDEBUG(__PRETTY_FUNCTION__ << "nothing to setup")
}



// void AlgorithmAddition::setTarget(PARTICLERANGE* range)
// {
//     vesDEBUG(__PRETTY_FUNCTION__)
//     target_range = enhance::make_observer<PARTICLERANGE>(range);
// }



void AlgorithmAddition::setParentSystem(System* sys)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    parent_system = enhance::make_observer<System>(sys);
}







void GrandCanonicalAddition::setup()
{
    distributor = std::make_unique<RandomDistributor>();
    distributor->setParameters(getParameters());
    distributor->check_for_aligned_box_setup();
    
    metropolis.setParameters(getParameters());
}



void GrandCanonicalAddition::apply()
{
    if(parent_system->getAlgorithm()->getCurrentStep() % 1000 == 0)
    {
        boost::progress_timer timer;
        auto a = tryRemove();
        auto b = tryInsert();
    }
}



bool GrandCanonicalAddition::tryRemove()
{

    return true;
}



bool GrandCanonicalAddition::tryInsert()
{
    vesLOG(__PRETTY_FUNCTION__);
    assert(distributor);
    vesLOG("calc potential energy before");
    auto energy_before = parent_system->potentialEnergy();
    vesLOG("  got " << energy_before);
    vesLOG("add particle");
    parent_system->addParticles(ParticleFactory<ParticleMobile>(1));

    auto new_particle = parent_system->getParticles().rbegin()->get();
    // new_particle->setParameters(getParameters());
    vesLOG("set rendom coords");
    new_particle->setCoords(distributor->randomCoords());
    vesLOG("set rendom orientation");
    new_particle->setOrientation(distributor->randomOrientation());

    while(distributor->conflicting_placement(&parent_system->getParticles(), *parent_system->getParticles().rbegin()))
    {
        vesLOG("placement confilct: set new coords and orientation");
        new_particle->setCoords(distributor->randomCoords());
        new_particle->setOrientation(distributor->randomOrientation());
    }

    vesLOG("calc potential energy after");
    auto energy_after = parent_system->potentialEnergy();
    vesLOG("  got " << energy_after);

    vesLOG("metropolis");
    if(metropolis.isValid(energy_after-energy_before))
    {
        vesLOG("added particle with d_energy " << energy_after-energy_before);
        parent_system->getAlgorithm()->getCells().deployParticle(*parent_system->getParticles().rbegin());
        vesLOG("deployed particle")
        return true;
    }
    else
    {
        vesLOG("no particle added with d_energy " << energy_after-energy_before);
        parent_system->getParticles().pop_back();
        return false;
    }


}