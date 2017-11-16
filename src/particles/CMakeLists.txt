set(particles_source_files
    particle_base.cpp
    particle_base.hpp
    particle_frame.hpp
    particle_mobile.cpp
    particle_mobile.hpp
)

find_package(Eigen3 REQUIRED NO_MODULE)
add_library(particles ${particles_source_files})

target_link_libraries(particles enhance Eigen3::Eigen)