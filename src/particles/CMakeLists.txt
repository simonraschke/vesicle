message(${CMAKE_CURRENT_SOURCE_DIR})

set(particles_source_files
    particle_base.cpp
    particle_mobile.cpp

    particle_base.hpp
    particle_mobile.hpp
    particle_frame.hpp

    particle_factory.hpp
)

find_package(Eigen3 REQUIRED NO_MODULE)
add_library(particles STATIC ${particles_source_files})

target_link_libraries(particles enhance Eigen3::Eigen)
target_include_directories(particles PUBLIC ${PROJECT_SOURCE_DIR})