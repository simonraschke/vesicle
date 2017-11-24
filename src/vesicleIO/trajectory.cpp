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



void TrajectoryWriter::setSkip(unsigned int s)
{
    skip = s;
}



bool TrajectoryWriter::isSkip()
{
    return skip_counter < skip;
}



void TrajectoryWriter::setAnisotropic(bool b)
{
    anisotropic = b;
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


void TrajectoryWriterGro::write(const HistoryStorage& history)
{
    if(likely(isSkip()))
    {
        ++skip_counter;
        return;
    }
    else
    {
        skip_counter = 0;
    }
    
    assert(FILE.is_open());
    assert(filename);
    assert(target_range);

    FILE << "t=" << history.getTime().back() << '\n';

    if(anisotropic)
        FILE << target_range->size()*2 << '\n';
    else
        FILE << target_range->size() << '\n';

    unsigned long atom = 0;
    for(unsigned long residue = 0; residue < target_range->size(); ++residue)
    {
        const auto& target = target_range->operator[](residue);
        const cartesian& coords = scaleDown(target->coords());
        const cartesian velocity = cartesian::Zero();

        FILE << std::setw(5) <<  residue+1;
        FILE << std::setw(5) <<  target->name();
        FILE << std::setw(5) <<  "A";
        FILE << std::setw(5) <<  atom+1;
        FILE << std::setprecision(3);
        FILE << std::setw(8) <<  coords(0);
        FILE << std::setw(8) <<  coords(1);
        FILE << std::setw(8) <<  coords(2);
        FILE << std::setprecision(4);
        FILE << std::setw(8) <<  velocity(0);
        FILE << std::setw(8) <<  velocity(1);
        FILE << std::setw(8) <<  velocity(2);
        FILE << '\n';
        ++atom;

        if(anisotropic)
        {
            const cartesian orientation = target->orientation().normalized();
            const cartesian& circularVelocity = target->circularVelocity();
            FILE << std::setw(5) <<  residue+1;
            FILE << std::setw(5) <<  target->name();
            FILE << std::setw(5) <<  "B";
            FILE << std::setw(5) <<  atom+1;
            FILE << std::setprecision(3);
            FILE << std::setw(8) <<  coords(0) + orientation(0)*getParameters().kappa/2.f;
            FILE << std::setw(8) <<  coords(1) + orientation(1)*getParameters().kappa/2.f;
            FILE << std::setw(8) <<  coords(2) + orientation(2)*getParameters().kappa/2.f;
            FILE << std::setprecision(4);
            FILE << std::setw(8) <<  circularVelocity(0);
            FILE << std::setw(8) <<  circularVelocity(1);
            FILE << std::setw(8) <<  circularVelocity(2);
            FILE << '\n';
            ++atom;
        }
    }

    FILE << getLengthX() << ' ' << getLengthY() << ' ' << getLengthZ();
    FILE << '\n';
}