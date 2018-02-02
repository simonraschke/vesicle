#include "gro_reader.hpp"



TrajectoryReaderGro::TrajectoryReaderGro()
    : TrajectoryReader()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}



void TrajectoryReaderGro::readAllFrames()
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(!FILE.is_open())
        throw std::logic_error("No open filestream to read from");
    else if (FILE.is_open())
    {
        Frame::first_type frame_counter = 0;
        Frame::second_type frame {};
        
        bool first_frame_header = true;

        for( std::string line; std::getline(FILE, line); /**/ )
        {
            if(boost::algorithm::contains(line, "FRAMEBEGIN"))
            {
                if( !first_frame_header )
                {
                    frames.insert(std::make_pair(frame_counter,frame));
                    ++frame_counter;
                    frame.clear();
                }
                first_frame_header = false;
            }
            frame.push_back(line);
        }

        // do not forget to add the last frame
        frames.insert(std::make_pair(frame_counter,frame));

        // std::move(std::begin(frame), std::end(frame), std::back_inserter(last_frame));
        vesLOG("last frame is frame " << frames.rbegin()->first )
    }
}



TrajectoryReaderGro::Frame TrajectoryReaderGro::getFrame(long long i) const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(frames.empty()) throw std::runtime_error("there is no frame to return");

    if(i < 0)
    {
        auto rit = std::crbegin(frames);
        std::advance(rit, std::abs(i)-1);
        return *rit;
    }
    else
    {
        auto it = std::cbegin(frames);
        std::advance(it, i);
        return *it;
    }
}



TrajectoryReaderGro::FrameMap TrajectoryReaderGro::getMatches(std::regex reg) const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(frames.empty()) throw std::runtime_error("there is no frame to return");

    FrameMap selection;

    if(std::regex_match("-1", reg))
    {
        vesLOG("found frame: " << frames.rbegin()->first)
        selection.insert(*frames.rbegin());
    }
    else
        for(const Frame& frame : frames)
        {
            // vesLOG("frame.first " << frame.first)
            if(std::regex_match(std::to_string(frame.first), reg))
            {
                vesLOG("found frame: " << frame.first)
                selection.insert(frame);
            }
        }

    if(selection.empty()) throw std::runtime_error("could not find any matching pattern while reading input frames");
    return selection;
}



const TrajectoryReaderGro::FrameMap& TrajectoryReaderGro::getFrames() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(frames.empty()) throw std::runtime_error("there is no frame to return");
    return frames;
}



std::size_t TrajectoryReaderGro::numParticles() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(frames.empty()) throw std::logic_error("no frame was read");
    std::string number = getFrame(0).second[1];
    return std::stoll( number ) / ( isAnisotropic() ? 2 : 1 );
}



bool TrajectoryReaderGro::isAnisotropic() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    if(frames.empty()) throw std::logic_error("no frame was read");
    std::string number_a { *( std::begin(getFrame(0).second[2]) + 14 ) };
    std::string number_b { *( std::begin(getFrame(0).second[3]) + 14 ) };
    try
    {
        return std::string( number_a ) != std::string( number_b );
    }
    catch (const std::exception& e)
    {
        vesCRITICAL("string comparison failed   " << number_a << " != " << number_b)
        throw;
    }
}



std::map<std::string,std::string> TrajectoryReaderGro::particleLineTokens(std::string line)
{
    vesDEBUG(__PRETTY_FUNCTION__)
    std::map<std::string,std::string> map;
    if(line.size() == 44)
    {
        line.append(std::string("  0.0000"));
        line.append(std::string("  0.0000"));
        line.append(std::string("  0.0000"));
    }

    map["resnum"] = line.substr(0,5);
    map["resname"] = line.substr(5,5);
    map["atomname"] = line.substr(10,5);
    map["atomnum"] = line.substr(15,5);
    map["pos x"] = line.substr(20,8);
    map["pos y"] = line.substr(28,8);
    map["pos z"] = line.substr(36,8);
    map["vel x"] = line.substr(44,8);
    map["vel y"] = line.substr(52,8);
    map["vel z"] = line.substr(60,8);

    for( auto& pair : map )
    #ifdef BOOST_VERSION
        boost::algorithm::erase_all(pair.second, " ");
    #else
        pair.second.erase(std::remove_if(std::begin(pair.second), std::end(pair.second), ::isspace), std::end(pair.second));
    #endif
    
    vesDEBUG("resnum " << map["resnum"] << " resname " << map["resname"] << " atomname " << map["atomname"] << " atomnum " << map["atomnum"])
    vesDEBUG("pos x " << map["pos x"] << " pos y " << map["pos y"] << "  pos z " << map["pos z"] << " vel x " << map["vel x"] << " vel y " << map["vel y"] << "  vel z " << map["vel z"])

    return map;
}