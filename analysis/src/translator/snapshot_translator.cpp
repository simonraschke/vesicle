#include "snapshot_translator.hpp"



void AnisotropicSnapshotTranslatorGro::operator()(TrajectoryReader::Frame snapshot) 
{
    particles.clear();
    ParticleIDGenerator().reset();
    num_particles = 0;
    time_elapsed = 0;
    x = y = z = 0;

    // const auto number_of_snapshot = snapshot.first;
    const auto lines = snapshot.second;

    assert(boost::algorithm::contains(lines.front(), "FRAMEBEGIN"));
    if(!boost::algorithm::contains(lines.front(), "FRAMEBEGIN"))
        throw std::runtime_error("first line of snapshot missing keyword");

    try
    {
        time_elapsed = boost::lexical_cast<float>(enhance::splitAtDelimiter(lines[0],"=").at(1));
        vesLOG("snapshot of time " << time_elapsed)
    }
    catch (const boost::bad_lexical_cast& e)
    {
        vesWARNING("lexical cast of " + enhance::splitAtDelimiter(lines[0],"=").at(1) + " resulted in")
        vesCRITICAL("" << e.what())
    }

    try
    {
        num_particles = boost::lexical_cast<std::size_t>(lines[1])/2;
        vesLOG("number of particles " << num_particles)
        if(particles.capacity() < num_particles)
            particles.reserve(num_particles);
    }
    catch (const boost::bad_lexical_cast& e)
    {
        vesWARNING("lexical cast of " + lines[1]+ " resulted in")
        vesCRITICAL("" << e.what())
    }
    
    enhance::IncrementalNumberGenerator<AnisotropicSnapshotTranslatorGro> ID_maker;
    for(std::size_t i = 0; i < num_particles; ++i)
    {
        const std::string line1 = lines[i*2+2];
        const std::string line2 = lines[i*2+3];
        TrajectoryReaderGro reader;
        auto tokens1 = reader.particleLineTokens(line1);
        auto tokens2 = reader.particleLineTokens(line2);

        if ( boost::algorithm::contains(tokens1.at("resname"), "MOBIL") && boost::algorithm::contains(tokens2.at("resname"), "MOBIL") )
        {
            ParticleFactory<ParticleMobile> factory(1);
            particles.emplace_back(factory.createParticle());
        }
        else if ( boost::algorithm::contains(tokens1.at("resname"), "FRAME") && boost::algorithm::contains(tokens2.at("resname"), "FRAME") )
        {
            ParticleFactory<ParticleFrame> factory(1);
            particles.emplace_back(factory.createParticle());
        }
        else
            throw std::logic_error(__PRETTY_FUNCTION__ + enhance::toStringViaStream(" resname are ", tokens1.at("resname")," and ", tokens2.at("resname")));
        
        particles.back()->setCoords(AnisotropicCoordsTranslatorGro()(line1,line2));
        particles.back()->setOrientation(AnisotropicOrientationTranslatorGro()(line1,line2));
    }
    
    const auto box_values = enhance::splitAtDelimiter(lines.back()," ");
    assert(box_values.size() == 3);
    x = boost::lexical_cast<float>(box_values[0]);
    y = boost::lexical_cast<float>(box_values[1]);
    z = boost::lexical_cast<float>(box_values[2]); 
}