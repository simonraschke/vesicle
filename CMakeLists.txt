project(vesicle)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_CXX_STANDARD 14)

# specific flags
set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g3 -ggdb -D_DEBUG -DDEBUG          -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -g0       -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

# general flags
add_compile_options(-Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++)
add_compile_options(-ftree-vectorize)
add_compile_options(-march=nocona -mtune=generic)

# Building the libraries
add_subdirectory(enhance)
add_subdirectory(src)

# Build executable from main
set(vesicle_source_files src/vesicle.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")
add_executable(vesicle ${vesicle_source_files})

# Link
set( CMAKE_EXE_LINKER_FLAGS "-ltbb -pthread -ltbbmalloc -lboost_system -lboost_filesystem -lboost_program_options" )
target_link_libraries( vesicle systems IO interactions ${CMAKE_EXE_LINKER_FLAGS} )

install( TARGETS vesicle RUNTIME DESTINATION bin )

enable_testing()
add_subdirectory(tests)