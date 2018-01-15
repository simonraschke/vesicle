#include "gro_reader.hpp"



TrajectoryReaderGro::TrajectoryReaderGro()
    : TrajectoryReader()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}



void TrajectoryReaderGro::readAllFrames()
{
    if(!FILE.is_open())
        throw std::logic_error("No open filestream to read from");
    else if (FILE.is_open())
    {
        std::deque<std::string> frame {};
        std::size_t frame_counter = 0;
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



std::map<TrajectoryReaderGro::Frame::first_type,TrajectoryReaderGro::Frame::second_type> TrajectoryReaderGro::getMatches(std::regex reg) const
{
    if(frames.empty()) throw std::runtime_error("there is no frame to return");

    std::map<TrajectoryReaderGro::Frame::first_type,TrajectoryReaderGro::Frame::second_type> selection;

    for(const Frame& frame : frames)
    {
        vesLOG("frame.first " << frame.first)
        if(std::regex_match(std::to_string(frame.first), reg))
        {
            vesLOG("found: " << frame.first)
            selection.insert(frame);
        }
    }

    if(selection.empty()) throw std::runtime_error("could not find any matching pattern while reading input frames");
    return selection;
}



const std::map<TrajectoryReaderGro::Frame::first_type,TrajectoryReaderGro::Frame::second_type>& TrajectoryReaderGro::getFrames() const
{
    if(frames.empty()) throw std::runtime_error("there is no frame to return");
    return frames;
}