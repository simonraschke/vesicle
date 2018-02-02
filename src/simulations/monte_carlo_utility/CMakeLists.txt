message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Eigen3 REQUIRED NO_MODULE)

file(GLOB monte_carlo_utility_source_files
    *.cpp
    *.hpp
)

add_library(mcutil STATIC ${monte_carlo_utility_source_files})

target_link_libraries(mcutil enhance particles Eigen3::Eigen "-ltbb -pthread -ltbbmalloc")
target_include_directories(mcutil PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(mcutil PUBLIC ${PROJECT_SOURCE_DIR}/src)