#pragma once

#include "definitions.hpp"
#include "parameters.hpp"
#include "enhance/observer_ptr.hpp"
#include "systems/box.hpp"
#include <boost/filesystem.hpp>
#include <iomanip>
#include <string>



struct TrajectoryWriter
    : public virtual ParameterDependentComponent
    , public virtual Box<PERIODIC::ON>
{
    typedef Particle::cartesian cartesian;
    typedef boost::filesystem::path PATH;
    typedef boost::filesystem::ifstream IFSTREAM;
    typedef boost::filesystem::ofstream OFSTREAM;

    // virtual void read(PATH) = 0;
    virtual void write() = 0;
    virtual void setFilename(std::string) = 0;

    void setTarget(PARTICLERANGE*);

    virtual ~TrajectoryWriter();

protected:
    using Box<PERIODIC::ON>::scaleDown;
    using Box<PERIODIC::ON>::getLengthX;
    using Box<PERIODIC::ON>::getLengthY;
    using Box<PERIODIC::ON>::getLengthZ;

    // virtual std::string format(const cartesian&) = 0;

    TrajectoryWriter();

    PATH working_dir;
    OFSTREAM FILE {};
    std::unique_ptr<std::string> filename {nullptr};
    enhance::observer_ptr<PARTICLERANGE> target_range {nullptr};
};



struct TrajectoryWriterGro
    : public TrajectoryWriter
{
    TrajectoryWriterGro();

    virtual void setFilename(std::string) override;
    virtual void write() override;   
};
