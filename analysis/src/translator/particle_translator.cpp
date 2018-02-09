#include "particle_translator.hpp"



Particle::cartesian AnisotropicCoordsTranslatorGro::operator()(std::string line1, std::string line2) const
{
    TrajectoryReaderGro reader;
    auto tokens1 = reader.particleLineTokens(line1);
    auto tokens2 = reader.particleLineTokens(line2);

    if (tokens1.at("atomname") == "B" && tokens2.at("atomname") == "A")
        ;
    else if (tokens1.at("atomname") == "A" && tokens2.at("atomname") == "B")
        std::swap(tokens1,tokens2);
    else
        throw std::logic_error(__PRETTY_FUNCTION__ + enhance::toStringViaStream(" atomnames are ", tokens1.at("atomname")," and ", tokens2.at("atomname")));

    Particle::cartesian position;
    {
        position(0) = (boost::lexical_cast<float>(tokens1.at("pos x")) + boost::lexical_cast<float>(tokens2.at("pos x"))) /2;
        position(1) = (boost::lexical_cast<float>(tokens1.at("pos y")) + boost::lexical_cast<float>(tokens2.at("pos y"))) /2;
        position(2) = (boost::lexical_cast<float>(tokens1.at("pos z")) + boost::lexical_cast<float>(tokens2.at("pos z"))) /2;
    }

    return position;
}



Particle::cartesian AnisotropicOrientationTranslatorGro::operator()(std::string line1, std::string line2) const
{
    TrajectoryReaderGro reader;
    auto tokens1 = reader.particleLineTokens(line1);
    auto tokens2 = reader.particleLineTokens(line2);

    if (tokens1.at("atomname") == "B" && tokens2.at("atomname") == "A")
        ;
    else if (tokens1.at("atomname") == "A" && tokens2.at("atomname") == "B")
        std::swap(tokens1,tokens2);
    else
        throw std::logic_error(__PRETTY_FUNCTION__ + enhance::toStringViaStream(" atomnames are ", tokens1.at("atomname")," and ", tokens2.at("atomname")));

    Particle::cartesian orientation;
    {
        Particle::cartesian position;
        position(0) = (boost::lexical_cast<float>(tokens1.at("pos x")) + boost::lexical_cast<float>(tokens2.at("pos x"))) /2;
        position(1) = (boost::lexical_cast<float>(tokens1.at("pos y")) + boost::lexical_cast<float>(tokens2.at("pos y"))) /2;
        position(2) = (boost::lexical_cast<float>(tokens1.at("pos z")) + boost::lexical_cast<float>(tokens2.at("pos z"))) /2;
        orientation(0) = boost::lexical_cast<float>(tokens1.at("pos x")) - position(0);
        orientation(1) = boost::lexical_cast<float>(tokens1.at("pos y")) - position(1);
        orientation(2) = boost::lexical_cast<float>(tokens1.at("pos z")) - position(2);
    }

    return orientation;
}



Particle* AnisotropicParticleTranslatorGro::operator()(std::string line1, std::string line2) const
{
    TrajectoryReaderGro reader;
    auto tokens1 = reader.particleLineTokens(line1);
    auto tokens2 = reader.particleLineTokens(line2);

    if ( boost::algorithm::contains(tokens1.at("resname"), "MOBIL") && boost::algorithm::contains(tokens2.at("resname"), "MOBIL") )
    {
        ParticleMobile* particle = new ParticleMobile();
        particle->setCoords(AnisotropicCoordsTranslatorGro()(line1,line2));
        particle->setOrientation(AnisotropicOrientationTranslatorGro()(line1,line2));
        return particle;
    }
    else
        throw std::logic_error(__PRETTY_FUNCTION__ + enhance::toStringViaStream(" resname are ", tokens1.at("resname")," and ", tokens2.at("resname")));
}