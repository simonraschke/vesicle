#include "gro_reader.hpp"



TrajectoryReaderGro::TrajectoryReaderGro()
    : TrajectoryReader()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}



void TrajectoryReaderGro::readLastFrame()
{
    if(!FILE.is_open())
        throw std::logic_error("No open filestream to read from");
    else if (FILE.is_open())
    {
        std::deque<std::string> frame;
        std::string line;
        bool is_last_frame = false;
        std::getline(FILE, line);

        while( !is_last_frame )
        {
            frame.clear();
            bool frame_complete = false;
            frame.push_back(line);

            while( std::getline(FILE, line) && !frame_complete )
            {
                if( boost::algorithm::contains(frame.back(), "FRAMEBEGIN") )
                {
                    frame_complete = true;
                }
                else 
                {
                    frame.push_back(line);
                }
            }
            if( FILE.eof() )
            {
                is_last_frame = true;
            }
        }
        last_frame = frame;
    }
}


void TrajectoryReaderGro::setFilename(std::string name)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << name)
    if(!filename && FILE.is_open())
    {
        FILE.close();
    }

    filename = std::make_unique<std::string>(name+".gro");

    FILE.open(*filename);
}