message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Boost COMPONENTS program_options system filesystem REQUIRED)

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

find_package(HDF5 REQUIRED)

file(GLOB vesicle_source_files
    */*.cpp
    */*/*.cpp
)

add_library(vesicle_lib STATIC ${vesicle_source_files})

set_target_properties(vesicle_lib PROPERTIES CXX_STANDARD 17)
target_compile_definitions(vesicle_lib PRIVATE 
    USE_MATH_DEFINES 
    EIGEN_FAST_MATH=1 
    EIGEN_MAX_CPP_VER=17
    EIGEN_DONT_PARALLELIZE
)
target_compile_options(vesicle_lib PRIVATE -Wall -Wextra -Wpedantic -march=nocona -std=c++17)

target_link_libraries(vesicle_lib enhance tbb pthread tbbmalloc boost_system boost_filesystem boost_program_options hdf5)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/enhance)
target_include_directories(vesicle_lib PUBLIC ${PROJECT_SOURCE_DIR}/lib/h5xx)