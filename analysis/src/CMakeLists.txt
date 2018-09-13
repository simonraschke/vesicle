message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})
set(_hdf5_libs hdf5)

find_package(VTK REQUIRED)
vtk_module_config(VTK
    vtkCommonCore # vtkPoints
    vtkFiltersGeometry #vtkDataSetSurfaceFilter
    vtkFiltersModeling #vtkSphereSource
    vtkIOXML #vtkXMLPolyDataWriter
)

file(GLOB vesicle_analysis_lib_source_files
    */*.cpp
    */*.cxx
)

add_library(vesicle_analysis_lib STATIC ${vesicle_analysis_lib_source_files})
set_target_properties(vesicle_analysis_lib PROPERTIES CXX_STANDARD 17)
target_compile_options(vesicle_analysis_lib PRIVATE -Wall -Wextra -Wpedantic -Weffc++ -march=nocona  -std=c++17)

target_link_libraries(vesicle_analysis_lib vesicle_lib enhance tbb pthread tbbmalloc ${VTK_LIBRARIES} ${_hdf5_libs})
target_include_directories(vesicle_analysis_lib PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(vesicle_analysis_lib PUBLIC ${vesicle_SOURCE_DIR}/src)
target_include_directories(vesicle_analysis_lib PUBLIC ${vesicle_SOURCE_DIR}/analysis/src)
target_include_directories(vesicle_analysis_lib PUBLIC ${vesicle_SOURCE_DIR}/analysis/HighFive/include/highfive)
target_include_directories(vesicle_analysis_lib PUBLIC ${VTK_INCLUDE_DIRS})