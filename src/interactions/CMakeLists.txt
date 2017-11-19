message(${CMAKE_CURRENT_SOURCE_DIR})

set(interactions_source_files
    lennard_jones.cpp

    lennard_jones.hpp
    interaction_interface.hpp
)

add_library(interactions STATIC ${interactions_source_files})

target_link_libraries(interactions particles)
target_include_directories(interactions PUBLIC ${PROJECT_SOURCE_DIR}/src)