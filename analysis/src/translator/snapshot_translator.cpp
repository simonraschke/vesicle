#include "snapshot_translator.hpp"



void AnisotropicSnapshotTranslatorGro::operator()(TrajectoryReader::Frame snapshot) 
{
    const auto number_of_snapshot = snapshot.first;
    const auto lines = snapshot.second;
    assert(boost::algorithm::contains(lines.front(), "FRAMEBEGIN"));
    time_elapsed = boost::lexical_cast<float>(enhance::splitAtDelimiter(lines[0],"t=").at(1));
    num_particles = boost::lexical_cast<std::size_t>(lines[0])/2;
    particles.reserve(num_particles);
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
        // else if ( boost::algorithm::contains(tokens1.at("resname"), "FRAME") && boost::algorithm::contains(tokens2.at("resname"), "FRAME") )
        // {
        //     ParticleFactory<ParticleFrame> factory(1);
        //     particles.emplace_back(factory.createParticle());
        // }
        else
            throw std::logic_error(__func__ + enhance::toStringViaStream("resname are ", tokens1.at("resname")," and ", tokens2.at("resname")));
        
        particles.back()->setCoords(AnisotropicCoordsTranslatorGro()(line1,line2));
        particles.back()->setOrientation(AnisotropicOrientationTranslatorGro()(line1,line2));
    }
    
    const auto box_values = enhance::splitAtDelimiter(lines.back()," ");
    assert(box_values.size() == 3);
    x = boost::lexical_cast<float>(box_values[0]);
    y = boost::lexical_cast<float>(box_values[1]);
    z = boost::lexical_cast<float>(box_values[2]); 
}