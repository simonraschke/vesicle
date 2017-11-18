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

set(CMAKE_CXX_STANDARD 14)

# add_library(enhance SHARED enhance/CMakeLists.txt)
# set_target_properties(enhance PROPERTIES LINKER_LANGUAGE CXX)

# Linking the libraries
add_subdirectory(enhance)
add_subdirectory(src)

# Build executable from main
set(vesicle_source_files src/vesicle.cpp)
add_executable(vesicle ${vesicle_source_files})

# build into bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")

# Link
set( CMAKE_EXE_LINKER_FLAGS "-ltbb" )
target_link_libraries( vesicle enhance particles simulations systems ${CMAKE_EXE_LINKER_FLAGS} )

install( TARGETS vesicle RUNTIME DESTINATION ~/bin )