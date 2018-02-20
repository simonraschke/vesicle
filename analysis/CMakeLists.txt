message(${CMAKE_CURRENT_SOURCE_DIR})

cmake_minimum_required(VERSION 3.0)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

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

set_property(TARGET vesicle_analysis PROPERTY CXX_STANARD 17)

add_subdirectory(HighFive)
add_subdirectory(src)

# Link
target_link_libraries(vesicle_analysis parsers tbb pthread tbbmalloc)

install(TARGETS vesicle_analysis RUNTIME DESTINATION bin)

enable_testing()