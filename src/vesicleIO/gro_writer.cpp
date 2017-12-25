#include "gro_writer.hpp"




TrajectoryWriterGro::TrajectoryWriterGro()
    : TrajectoryWriter()
{
    vesDEBUG(__PRETTY_FUNCTION__)
}



void TrajectoryWriterGro::setAnisotropic(bool b)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << b)
    anisotropic = b;
    makeStartFileVMD();
}



void TrajectoryWriterGro::setFilename(std::string name)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << name)
    if(!filename && FILE.is_open())
    {
        FILE.close();
    }

    filename = std::make_unique<std::string>(name+".gro");

    FILE.open(*filename);
    FILE.setf(std::ios::fixed | std::ios::showpoint);

    makeStartFileVMD();
}


void TrajectoryWriterGro::write(const HistoryStorage& history)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  simulation time: " << history.getTime().back())
    ++skip_counter;
    if(skip_counter%getParameters().traj_skip!=0)
    {
        return;
    }
    else
    {
        skip_counter = 1;
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
        const cartesian& coords = scaleDownForVMD(target->coords());

        if(!anisotropic)
        {
            const cartesian velocity = target->velocity();
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
        }
        else
        {
            const cartesian orientation = target->orientation().normalized();
            const cartesian displacement = target->orientation() - target->orientationOld();
            
            FILE << std::setw(5) <<  residue+1;
            FILE << std::setw(5) <<  target->name();
            FILE << std::setw(5) <<  "A";
            FILE << std::setw(5) <<  atom+1;
            FILE << std::setprecision(3);
            FILE << std::setw(8) <<  coords(0) - orientation(0)*getParameters().kappa/2.f;
            FILE << std::setw(8) <<  coords(1) - orientation(1)*getParameters().kappa/2.f;
            FILE << std::setw(8) <<  coords(2) - orientation(2)*getParameters().kappa/2.f;
            FILE << std::setprecision(4);
            FILE << std::setw(8) <<  displacement(0);
            FILE << std::setw(8) <<  displacement(1);
            FILE << std::setw(8) <<  displacement(2);
            FILE << '\n';
            ++atom;

            FILE << std::setw(5) <<  residue+1;
            FILE << std::setw(5) <<  target->name();
            FILE << std::setw(5) <<  "B";
            FILE << std::setw(5) <<  atom+1;
            FILE << std::setprecision(3);
            FILE << std::setw(8) <<  coords(0) + orientation(0)*getParameters().kappa/2.f;
            FILE << std::setw(8) <<  coords(1) + orientation(1)*getParameters().kappa/2.f;
            FILE << std::setw(8) <<  coords(2) + orientation(2)*getParameters().kappa/2.f;
            FILE << std::setprecision(4);
            FILE << std::setw(8) <<  -displacement(0);
            FILE << std::setw(8) <<  -displacement(1);
            FILE << std::setw(8) <<  -displacement(2);
            FILE << '\n';
            ++atom;
        }
    }

    FILE << getLengthX() << ' ' << getLengthY() << ' ' << getLengthZ();
    FILE << '\n';
}



void TrajectoryWriterGro::makeStartFileVMD() const
{
    vesDEBUG(__PRETTY_FUNCTION__)
    // OFSTREAM STARTER;
    // STARTER.open(vesicle.vmd);
    // STARTER << "#!/bin/bash" << '\n';
    // STARTER << "vmd" << '\n';
    // STARTER.close();

    OFSTREAM VMD;
    VMD.open(".vmdrc");
    VMD << "mol load gro " << *filename << '\n';
    VMD << "light 0 on" << '\n';
    VMD << "light 1 on" << '\n';
    VMD << "light 2 on" << '\n';
    VMD << "light 3 on" << '\n';
    VMD << "display nearclip set 0" << '\n';
    
    VMD << "axes location off" << '\n';
    VMD << "stage location off" << '\n';
    
    VMD << "menu main on" << '\n';
    VMD << "menu graphics on" << '\n';
    
    VMD << "display resize 1080 1080" << '\n';
    VMD << "display reposition" << '\n';
    VMD << "display projection perspective" << '\n';
    VMD << "display rendermode GLSL" << '\n';
    VMD << "color Display Background white" << '\n';
    // VMD << "draw color black" << '\n';
    VMD << "mol coloring 7 2 ResName" << '\n';
    VMD << "mol modstyle 0 0 Licorice 3 15 15" << '\n';
    VMD << "mol modmaterial 0 0 AOChalky" << '\n';
    
    VMD << "# color definitions" << '\n';
    VMD << "color change rgb  0 0.07 0.20 0.48 ;# blue" << '\n';
    VMD << "color change rgb  1 0.70 0.20 0.10 ;# red" << '\n';
    VMD << "color change rgb  2 0.40 0.40 0.40 ;# gray" << '\n';
    VMD << "color change rgb  3 0.70 0.40 0.00 ;# orange" << '\n';
    VMD << "color change rgb  4 0.80 0.70 0.10 ;# yellow" << '\n';
    VMD << "color change rgb  7 0.13 0.47 0.04 ;# green" << '\n';
    VMD << "color change rgb  8 1.00 1.00 1.00 ;# white" << '\n';
    VMD << "color change rgb 10 0.10 0.70 0.80 ;# cyan" << '\n';
    VMD << "color change rgb 11 0.60 0.10 0.60 ;# purple" << '\n';
    VMD << "color change rgb 16 0.15 0.15 0.15 ;# black" << '\n';

    VMD << "after idle {" << '\n';
    VMD << "  pbc box -center origin -color black -width 1" << '\n';
    VMD << "  # set colors" << '\n';
    VMD << "  # create dummy molecule with one atom" << '\n';
    VMD << "  set mol [mol new atoms 1]" << '\n';
    VMD << "  set sel [atomselect $mol all]" << '\n';
    VMD << "  # add items to color categories" << '\n';
    VMD << "  $sel set name A" << '\n';
    VMD << "  $sel set type A" << '\n';
    VMD << "  $sel set name B" << '\n';
    VMD << "  $sel set type B" << '\n';
    VMD << "  $sel set name C" << '\n';
    VMD << "  $sel set type C" << '\n';
    VMD << "  # now we can define colors" << '\n';
    VMD << "  color Name A black" << '\n';
    VMD << "  color Type A black" << '\n';
    VMD << "  color Name B 23" << '\n';
    VMD << "  color Type B 23" << '\n';
    VMD << "  color Name C gray" << '\n';
    VMD << "  color Type C gray" << '\n';
    VMD << "  mol delete $mol" << '\n';

    if(anisotropic)
    {
        VMD << "  # clean up" << '\n';
        VMD << "  $sel delete" << '\n';
        // VMD << "  set mol [mol new atoms 1]" << '\n';
        // VMD << "  set sel [atomselect $mol all]" << '\n';
        VMD << "  for {set x 0} {$x < " << target_range->size()*2 <<"} {incr x} {" << '\n';
        VMD << "    set y [expr $x+1]" << '\n'; 
        VMD << "    set sel [atomselect top \"index $x $y\"]" << '\n';
        VMD << "    incr x 1" << '\n';
        VMD << "    set bonds [$sel getbonds]" << '\n';
        VMD << "    set ids [$sel get index]" << '\n';
        VMD << "    lassign $bonds atom1bonds atom2bonds" << '\n';
        VMD << "    lassign $ids atom1id atom2id" << '\n';
        VMD << "    lappend atom1bonds $atom2id" << '\n';
        VMD << "    lappend atom2bonds $atom1id" << '\n';
        VMD << "    $sel setbonds [list $atom1bonds $atom2bonds]"  << '\n';
        VMD << "    }" << '\n';
    }
    VMD << "}" << '\n';
    
    VMD.close();
}