cmake_minimum_required(VERSION 3.1)
project(vesicle)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_CXX_STANDARD 14)

# Building the libraries
add_subdirectory(enhance)
add_subdirectory(src)
add_subdirectory(analysis)



# specific flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

### vesicle
set(vesicle_source_files src/vesicle.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")
add_executable(vesicle ${vesicle_source_files})

target_compile_options(vesicle PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++ -march=nocona -mtune=generic)

target_compile_definitions(vesicle PUBLIC 
    USE_MATH_DEFINES 
    EIGEN_FAST_MATH=1 
    EIGEN_MAX_CPP_VER=14 
    EIGEN_DONT_PARALLELIZE
)

set_target_properties(vesicle PROPERTIES CXX_STANARD 14)

# Link
set(CMAKE_EXE_LINKER_FLAGS "-ltbb -pthread -ltbbmalloc")
target_link_libraries(vesicle systems IO interactions)

install( TARGETS vesicle RUNTIME DESTINATION bin )

enable_testing()
add_subdirectory(tests)
