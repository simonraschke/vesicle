set(particles_source_files
    particle_base.hpp
    particle_frame.hpp
    particle_mobile.hpp
)

add_library(particles ${particles_source_files})

target_link_libraries(particles enhance)