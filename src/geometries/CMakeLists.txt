message(${CMAKE_CURRENT_SOURCE_DIR})

set(geometries_source_files
    grid.cpp

    geometry.hpp
    grid.hpp
)

find_package(Eigen3 REQUIRED NO_MODULE)
add_library(geometries STATIC ${geometries_source_files})

target_link_libraries(geometries Eigen3::Eigen)
target_include_directories(geometries PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(geometries PUBLIC ${PROJECT_SOURCE_DIR}/src)