#include "h5_writer.hpp"



TrajectoryWriterH5::TrajectoryWriterH5()
    : TrajectoryWriter()
{
    vesDEBUG(__PRETTY_FUNCTION__);

}



void TrajectoryWriterH5::setAnisotropic(bool b)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << b)
    anisotropic = b;
    makeStartFileVMD();
}



void TrajectoryWriterH5::setPath(PATH path)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << path)
    if(!file_path && FILE.is_open())
    {
        FILE.close();
    }

    file_path = std::make_unique<PATH>(boost::filesystem::system_complete(working_dir/path));

    if(boost::filesystem::exists(*file_path))
    {
        // splitting the filepath
        auto string_parts = enhance::splitAtDelimiter(file_path->string(), ".");

        std::string first_part = std::accumulate(string_parts.begin(), std::next(string_parts.rbegin()).base(), std::string(""), [](auto i, auto j){return i+j+".";});
        {
            std::string name_appendix = "_old";
            first_part.insert(std::next(first_part.rbegin()).base(), std::begin(name_appendix), std::end(name_appendix));
        }
        std::string filetype = *string_parts.rbegin();

        // making "trajectory_old.gro" from "trajectory.gro"
        PATH destination = boost::filesystem::system_complete(first_part+filetype);
        vesLOG("trajectory file " << file_path->string() << " already exists. will backup to " << destination.string())

        // and backup the old trajectory
        boost::filesystem::copy_file(*file_path, destination, boost::filesystem::copy_option::overwrite_if_exists);
    }

    if(GLOBAL::getInstance().mode == GLOBAL::NEWRUN)
    {
        vesLOG("trunc open " << file_path->string())
        h5file.open(file_path->string(), h5xx::file::trunc);
    }
    else
    {
        vesLOG("app open" << file_path->string())
        h5file.open(file_path->string(), h5xx::file::out);
    }
}



void TrajectoryWriterH5::write(const float& time_elapsed, bool FORCE)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  simulation time: " << time_elapsed);

    if(!FORCE)
    {
        ++skip_counter;
        if(skip_counter%getParameters().out_traj_skip!=0)
            return;
        else
            skip_counter = 0;
    }
    
    assert(h5file.valid());
    assert(file_path);
    assert(target_range);

    const std::string group_name("snapshot"+std::to_string(static_cast<std::size_t>(time_elapsed)));
    h5xx::group group(h5file, group_name);

    {
        h5xx::write_attribute(group, "general.algorithm", getParameters().algorithm);
        h5xx::write_attribute(group, "general.acceptance", getParameters().acceptance);
        h5xx::write_attribute(group, "general.interaction", getParameters().interaction);
        h5xx::write_attribute(group, "general.thermostat", getParameters().thermostat);
        h5xx::write_attribute(group, "general.grand_canonical", getParameters().grand_canonical);
        h5xx::write_attribute(group, "system.frame_guides_grid_edge", getParameters().frame_guides_grid_edge);
        h5xx::write_attribute(group, "system.guiding_elements_each", getParameters().guiding_elements_each);
        h5xx::write_attribute(group, "system.guiding_elements_plane", getParameters().guiding_elements_plane);
        h5xx::write_attribute(group, "system.plane_edge", getParameters().plane_edge);
        h5xx::write_attribute(group, "system.osmotic", getParameters().osmotic);
        h5xx::write_attribute(group, "system.osmotic_density_inside", getParameters().osmotic_density_inside);
        h5xx::write_attribute(group, "system.num_mobile_particles", std::count_if(target_range->begin(), target_range->end(), [](const auto& p){ return p->getType() == PARTICLETYPE::MOBILE;} ));
        h5xx::write_attribute(group, "system.num_frame_particles", std::count_if(target_range->begin(), target_range->end(), [](const auto& p){ return p->getType() == PARTICLETYPE::FRAME;} ));
        h5xx::write_attribute(group, "system.num_osmotic_particles", std::count_if(target_range->begin(), target_range->end(), [](const auto& p){ return p->getType() == PARTICLETYPE::OSMOTIC;} ));
        h5xx::write_attribute(group, "system.num_all_particles", target_range->size());
        h5xx::write_attribute(group, "system.box.x", getParameters().x);
        h5xx::write_attribute(group, "system.box.y", getParameters().y);
        h5xx::write_attribute(group, "system.box.z", getParameters().z);
        h5xx::write_attribute(group, "system.actual_time", time_elapsed);
        h5xx::write_attribute(group, "system.timestep", getParameters().dt);
        h5xx::write_attribute(group, "system.time_max", getParameters().time_max);
        h5xx::write_attribute(group, "system.temperature", getParameters().temperature);
        h5xx::write_attribute(group, "system.kappa", getParameters().kappa);
        h5xx::write_attribute(group, "system.gamma", getParameters().gamma);
        h5xx::write_attribute(group, "system.ljepsilon", getParameters().LJepsilon);
        h5xx::write_attribute(group, "system.ljsigma", getParameters().LJsigma);
        h5xx::write_attribute(group, "system.sw_position_min", getParameters().sw_position_min);
        h5xx::write_attribute(group, "system.sw_position_max", getParameters().sw_position_max);
        h5xx::write_attribute(group, "system.sw_position_target", getParameters().sw_position_target);
        h5xx::write_attribute(group, "system.sw_orientation_min", getParameters().sw_orientation_min);
        h5xx::write_attribute(group, "system.sw_orientation_max", getParameters().sw_orientation_max);
        h5xx::write_attribute(group, "system.sw_orientation_target", getParameters().sw_orientation_target);
        h5xx::write_attribute(group, "system.cell_min_edge", getParameters().cell_min_edge);
        h5xx::write_attribute(group, "system.max_cells_dim", getParameters().max_cells_dim);
        h5xx::write_attribute(group, "output.out_traj", getParameters().out_traj);
        h5xx::write_attribute(group, "output.skip", getParameters().out_traj_skip);
    }
    
    {
        const std::string dataset_name("position");
        auto positions = getPositions();
        h5xx::create_dataset(group, dataset_name, positions);
        h5xx::write_dataset(group, dataset_name, positions);
    }
    
    {
        const std::string dataset_name("orientation");
        auto orientations = getPositions();
        h5xx::create_dataset(group, dataset_name, orientations);
        h5xx::write_dataset(group, dataset_name, orientations);
    }
    
    {
        const std::string dataset_name("type");
        auto orientations = getResidueTypes();
        h5xx::create_dataset(group, dataset_name, orientations);
        h5xx::write_dataset(group, dataset_name, orientations);
    }
}



TrajectoryWriterH5::array2d_t TrajectoryWriterH5::getPositions() const
{
    typedef array2d_t::index index;

    array2d_t array(boost::extents[target_range->size()][3]);

    for(index residue = 0; residue < static_cast<index>(target_range->size()); ++residue)
    {
        const auto& target = target_range->operator[](residue);
        const cartesian& coords = scaleDown( target->coords());

        array[residue][static_cast<index>(0)] = coords(0);
        array[residue][static_cast<index>(1)] = coords(1);
        array[residue][static_cast<index>(2)] = coords(2);
    }
    
    return array;
}



TrajectoryWriterH5::array2d_t TrajectoryWriterH5::getOrientations() const
{
    typedef array2d_t::index index;

    array2d_t array(boost::extents[target_range->size()][3]);

    for(index residue = 0; residue < static_cast<index>(target_range->size()); ++residue)
    {
        const auto& target = target_range->operator[](residue);
        const cartesian orientation = target->getType() == OSMOTIC ? cartesian(0,0,0) : target->orientation().normalized();

        array[residue][static_cast<index>(0)] = orientation(0);
        array[residue][static_cast<index>(1)] = orientation(1);
        array[residue][static_cast<index>(2)] = orientation(2);
    }
    
    return array;
}



TrajectoryWriterH5::array1d_t TrajectoryWriterH5::getResidueTypes() const
{
    typedef array1d_t::index index;

    array1d_t array(boost::extents[target_range->size()]);

    for(index residue = 0; residue < static_cast<index>(target_range->size()); ++residue)
    {
        array[residue] = target_range->operator[](residue)->getType();
    }
    
    return array;
}







void TrajectoryWriterH5::makeStartFileVMD() const
{
    ;
}