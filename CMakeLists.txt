project(vesicle)
cmake_minimum_required(VERSION 3.0)

# add_library(enhance SHARED enhance/CMakeLists.txt)
# set_target_properties(enhance PROPERTIES LINKER_LANGUAGE CXX)

# Linking the libraries
add_subdirectory(src)
add_subdirectory(enhance)

# Where the header files are
include_directories(${CMAKE_SOURCE_DIR}/enhance)

# Build executable from main
set(vesicle_source_files src/vesicle.cpp)
add_executable(vesicle ${vesicle_source_files})

# build into bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "~/bin/")

# Link
set( CMAKE_EXE_LINKER_FLAGS "" )
target_link_libraries( vesicle particles ${CMAKE_EXE_LINKER_FLAGS} )