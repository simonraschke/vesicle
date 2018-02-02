message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Eigen3 REQUIRED NO_MODULE)

file(GLOB particles_source_files
    *.cpp
    *.hpp
)

add_library(particles STATIC ${particles_source_files})

target_link_libraries(particles enhance geometries IO Eigen3::Eigen)
target_include_directories(particles PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(particles PUBLIC ${PROJECT_SOURCE_DIR}/src)