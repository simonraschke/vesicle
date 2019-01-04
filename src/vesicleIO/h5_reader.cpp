#include "h5_reader.hpp"



TrajectoryReaderH5::TrajectoryReaderH5()
    : TrajectoryReader()
{
    vesDEBUG(__PRETTY_FUNCTION__);
}



void TrajectoryReaderH5::setPath(PATH path)
{
    // vesDEBUG(__PRETTY_FUNCTION__<< "  " << path)
    // if(!file_path && FILE.is_open())
    // {
    //     FILE.close();
    // }

    // file_path = std::make_unique<PATH>(boost::filesystem::system_complete(working_dir/path));

    // if(boost::filesystem::exists(*file_path))
    // {
    //     vesLOG("read open " << file_path->string())
    //     h5file.open(getParameters().in_traj_path.string(), h5xx::file::in);
    // }
}




void TrajectoryReaderH5::readAllFrames(bool direct_match_check)
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    // if(!h5file.valid())
    //     throw std::logic_error("No hdf5 file to read from");
    // else if(h5file.valid())
    // {
    //     frame_counter = 0;
    //     frame_buffer.clear();

    //     for(auto name : getGroupNames())
    //         std::cout << name << std::endl;
        

    //     // do not forget to add the last frame
    //     // directly check for regex?
    //     if(direct_match_check)
    //     {
    //         // special case if only the last frame is wanted
    //         if(std::regex_match("-1", getParameters().in_frames))
    //         {
    //             frames.emplace_back(std::make_pair(frame_counter,frame_buffer));
    //         }
    //         // check for regex match if no special case
    //         else if(isRegexMatch(std::make_pair(frame_counter,frame_buffer), getParameters().in_frames))
    //         {
    //             frames.emplace_back(std::make_pair(frame_counter,frame_buffer));
    //         }
    //     }
    //     else
    //     {
    //         frames.emplace_back(std::make_pair(frame_counter,frame_buffer));
    //     }

    //     vesLOG("last frame is frame " << frames.rbegin()->first )
    // }
}




void TrajectoryReaderH5::readNextFrame(std::regex reg)
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    // if(!FILE.is_open())
    //     throw std::logic_error("No open filestream to read from");
    // else if (FILE.eof())
    // {
    //     vesLOG("reached EOF")
    // }
    // else if (FILE.is_open())
    // {
    //     FILE.seekg(FILE_pos);

    //     frames.clear();
    //     frame_buffer.clear();
        
    //     bool first_frame_header = true;
    //     bool found = false;

    //     // iterating the whole file
    //     for( std::string line; std::getline(FILE, line) || FILE.eof(); FILE_pos = FILE.tellg() )
    //     {
    //         // if FRAMEBEGIN is in line
    //         if(boost::algorithm::contains(line, "FRAMEBEGIN") || FILE.eof())
    //         {
    //             // and its not the first line
    //             if( !first_frame_header )
    //             {
    //                 // add the whole frame_buffer to container

    //                 // directly check for regex?
    //                 if(isRegexMatch(std::make_pair(frame_counter,frame_buffer),reg))
    //                 {
    //                     found = true;
    //                     frames.emplace_back(std::make_pair(frame_counter,frame_buffer));
    //                     vesLOG("appending frame " << frame_counter)
    //                 }

    //                 // increase frame counter
    //                 ++frame_counter;

    //                 // clear the actual frame_buffer, because if you encounter
    //                 // FRAMEBEGIN the next frame is reached
    //                 frame_buffer.clear();
    //             }
    //             //first occurance of FRAMEBEGIN was obviously met
    //             first_frame_header = false;
    //         }
    //         // append line to frame_buffer
    //         frame_buffer.push_back(line);
    //         if(found)
    //             break;
    //         else 
    //             FILE_pos = FILE.tellg();
    //     }

    //     if(!found)
    //     {
    //         vesWARNING("no frame for regex was found")
    //     }
    //     else 
    //     {
    //         vesLOG("last frame is frame " << frames.rbegin()->first )
    //     }
    // }
}
