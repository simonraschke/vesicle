SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "g++")
cmake_minimum_required(VERSION 3.1)
project(vesicle)

message(WARNING "BUILD TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_VERBOSE_MAKEFILE OFF)

# Building the libraries
add_subdirectory(enhance)
add_subdirectory(src)
add_subdirectory(analysis)

# specific flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g3 -ggdb -D_DEBUG -DDEBUG -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g1 -ggdb -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g0 -DNDEBUG -DEIGEN_NO_DEBUG -ffast-math")

### vesicle
set(basic_source_files src/vesicle.cpp)
add_executable(vesicle ${basic_source_files})

set_target_properties(vesicle PROPERTIES RUNTIME_OUTPUT_DIRECTORY "~/bin/")
set_target_properties(vesicle PROPERTIES CXX_STANDARD 17)
target_compile_options(vesicle PRIVATE -Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++ -march=nocona -std=c++1z)

# Link
target_link_libraries(vesicle vesicle_lib)

install( TARGETS vesicle RUNTIME DESTINATION bin )

enable_testing()
add_subdirectory(tests)
