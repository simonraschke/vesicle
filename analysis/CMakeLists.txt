message(${CMAKE_CURRENT_SOURCE_DIR})

cmake_minimum_required(VERSION 3.0)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG -frounding-math -fsignaling-nans -fsanitize=leak")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math -fsanitize=leak")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

set(vesicle_analysis_source_files src/vesicle_analysis.cpp)
add_executable(vesicle_analysis ${vesicle_analysis_source_files})

set_target_properties(vesicle_analysis PROPERTIES RUNTIME_OUTPUT_DIRECTORY "~/bin/")
set_target_properties(vesicle_analysis PROPERTIES CXX_STANDARD 17)
target_compile_options(vesicle_analysis PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++ -march=nocona -std=c++17)
target_compile_definitions(vesicle_analysis PUBLIC 
    USE_MATH_DEFINES
    EIGEN_FAST_MATH=1 
    EIGEN_MAX_CPP_VER=14
    EIGEN_DONT_PARALLELIZE
)

add_subdirectory(HighFive)
add_subdirectory(src)

# Link
target_link_libraries(vesicle_analysis vesicle_analysis_lib)

install(TARGETS vesicle_analysis RUNTIME DESTINATION bin)

enable_testing()