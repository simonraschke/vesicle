project(vesicle)
cmake_minimum_required(VERSION 3.0)

# specific flags
set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g3 -ggdb -D_DEBUG -DDEBUG          -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -g0       -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

# general flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=nocona -mtune=generic")

set(CMAKE_VERBOSE_MAKEFILE Off)
set(CMAKE_CXX_STANDARD 14)

# Building the libraries
add_subdirectory(enhance)
add_subdirectory(src)

# Build executable from main
set(vesicle_source_files src/vesicle.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")
add_executable(vesicle ${vesicle_source_files})

# Link
set( CMAKE_EXE_LINKER_FLAGS "-ltbb -lboost_system -lboost_filesystem -lboost_program_options" )
target_link_libraries( vesicle systems IO interactions ${CMAKE_EXE_LINKER_FLAGS} )

install( TARGETS vesicle RUNTIME DESTINATION bin )