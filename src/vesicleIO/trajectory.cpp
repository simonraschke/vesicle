#include "trajectory.hpp"



TrajectoryWriter::TrajectoryWriter()
    : working_dir(boost::filesystem::current_path())
{

}



TrajectoryWriter::~TrajectoryWriter()
{
    if(FILE.is_open())
    {
        FILE.close();
    }
}



void TrajectoryWriter::setTarget(PARTICLERANGE* range)
{
    target_range = enhance::make_observer<PARTICLERANGE>(range);
}



TrajectoryWriterGro::TrajectoryWriterGro()
    : TrajectoryWriter()
{
    
}




void TrajectoryWriterGro::setFilename(std::string name)
{
    if(!filename && FILE.is_open())
    {
        FILE.close();
    }

    filename = std::make_unique<std::string>(name+".gro");

    FILE.open(*filename);
    FILE.setf(std::ios::fixed | std::ios::showpoint);
}


void TrajectoryWriterGro::write()
{
    assert(FILE.is_open());
    assert(filename);
    assert(target_range);

    FILE << "PLACEHOLDER" << '\n';
    FILE << target_range->size() << '\n';

    for(unsigned long i = 0; i < target_range->size(); ++i)
    {
        const auto& target = target_range->operator[](i);
        const cartesian& coords = scaleDown(target->coords());
        const cartesian& velocity = target->velocity();

        FILE << std::setw(5) <<  i+1;
        FILE << std::setw(5) <<  target->name();
        FILE << std::setw(5) <<  "A";
        FILE << std::setw(5) <<  i+1;
        FILE << std::setprecision(3);
        FILE << std::setw(8) <<  coords(0);
        FILE << std::setw(8) <<  coords(1);
        FILE << std::setw(8) <<  coords(2);
        FILE << std::setprecision(4);
        FILE << std::setw(8) <<  velocity(0);
        FILE << std::setw(8) <<  velocity(1);
        FILE << std::setw(8) <<  velocity(2);
        FILE << '\n';
    }

    FILE << getLengthX() << ' ' << getLengthY() << ' ' << getLengthZ();
    FILE << '\n';
}