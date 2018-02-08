message(${CMAKE_CURRENT_SOURCE_DIR})

cmake_minimum_required(VERSION 3.0)

add_subdirectory(HighFive)
add_subdirectory(src)

### vesicle_analysis
set(vesicle_analysis_source_files src/vesicle_analysis.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")
add_executable(vesicle_analysis ${vesicle_analysis_source_files})

target_compile_options(vesicle_analysis PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3 -mavx -Weffc++ -march=nocona -mtune=generic)

target_compile_definitions(vesicle_analysis PUBLIC 
    USE_MATH_DEFINES 
    EIGEN_FAST_MATH=1 
    EIGEN_MAX_CPP_VER=14 
    EIGEN_DONT_PARALLELIZE
)

set_target_properties(vesicle_analysis PROPERTIES CXX_STANARD 14 )

target_include_directories(vesicle_analysis PRIVATE 
    ${CMAKE_CURRENT_LIST_DIR}/HighFive/include/highfive
    ${vesicle_SOURCE_DIR}/src
)

# Link
target_link_libraries(vesicle_analysis translator systems tbb pthread tbbmalloc)

install(TARGETS vesicle_analysis RUNTIME DESTINATION bin)

enable_testing()