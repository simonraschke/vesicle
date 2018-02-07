# set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
# set(CMAKE_DISABLE_SOURCE_CHANGES  ON)


cmake_minimum_required(VERSION 3.1)
project(vesicle)

set(CMAKE_VERBOSE_MAKEFILE OFF)

# Building the libraries
add_subdirectory(enhance)
add_subdirectory(src)

# add_subdirectory(analysis)


# specific flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

### vesicle
set(vesicle_source_files src/vesicle.cpp)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")
add_executable(vesicle ${vesicle_source_files})
target_compile_options(vesicle PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3 -mavx -Weffc++ -march=nocona -mtune=generic)
target_compile_definitions(vesicle PUBLIC USE_MATH_DEFINES EIGEN_FAST_MATH=1 EIGEN_MAX_CPP_VER=14 EIGEN_DONT_PARALLELIZE)
set_target_properties(vesicle PROPERTIES 
    CXX_STANARD 14 
    CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/"
)

# Link
set(CMAKE_EXE_LINKER_FLAGS "-ltbb -pthread -ltbbmalloc")
target_link_libraries(vesicle systems IO interactions)

install( TARGETS vesicle RUNTIME DESTINATION bin )

# include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
# ExternalProject_Add(HighFive SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/HighFive PREFIX HighFive)

enable_testing()
add_subdirectory(tests)
<<<<<<< Updated upstream
# add_subdirectory(HighFive)
=======
add_subdirectory(HighFive)
>>>>>>> Stashed changes
