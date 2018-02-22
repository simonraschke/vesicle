message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB particles_source_files
    *.cpp
    *.hpp
)

add_library(particles STATIC ${particles_source_files})

set_target_properties(particles PROPERTIES CXX_STANARD 14)

# target_link_libraries(particles enhance geometries IO Eigen::Eigen)
target_link_libraries(particles enhance geometries IO)
target_include_directories(particles PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(particles PUBLIC ${PROJECT_SOURCE_DIR}/src)