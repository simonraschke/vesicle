message(${CMAKE_CURRENT_SOURCE_DIR})

cmake_minimum_required(VERSION 3.0)


### vesicle_analysis
set(vesicle_analysis_source_files src/vesicle_analysis.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")
add_executable(vesicle_analysis ${vesicle_analysis_source_files})

target_compile_options(vesicle_analysis PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++ -march=nocona -mtune=generic)

target_compile_definitions(vesicle_analysis PUBLIC 
USE_MATH_DEFINES 
EIGEN_FAST_MATH=1 
EIGEN_MAX_CPP_VER=14 
EIGEN_DONT_PARALLELIZE
)

set_target_properties(vesicle_analysis PROPERTIES CXX_STANARD 14 )

add_subdirectory(HighFive)
add_subdirectory(src)

# Link
target_link_libraries(vesicle_analysis parsers tbb pthread tbbmalloc)

install(TARGETS vesicle_analysis RUNTIME DESTINATION bin)

enable_testing()