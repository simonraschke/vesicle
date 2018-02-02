message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Boost COMPONENTS program_options system filesystem REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})
set(_hdf5_libs hdf5 hdf5_cpp)

file(GLOB IO_source_files
    *.cpp
    *.hpp
)

add_library(IO STATIC ${IO_source_files})

target_link_libraries(IO enhance Eigen3::Eigen ${_hdf5_libs} "-lboost_system -lboost_filesystem -lboost_program_options")

target_include_directories(IO PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(IO PUBLIC ${PROJECT_SOURCE_DIR}/src)