message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Eigen3 REQUIRED NO_MODULE)

set(interactions_source_files
    lennard_jones.cpp
    angular_lennard_jones.cpp

    lennard_jones.hpp
    angular_lennard_jones.hpp
    interaction.hpp
)

add_library(interactions STATIC  ${interactions_source_files})

target_link_libraries(interactions particles Eigen3::Eigen)
target_include_directories(interactions PUBLIC ${PROJECT_SOURCE_DIR}/src)