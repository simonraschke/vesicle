message(${CMAKE_CURRENT_SOURCE_DIR})

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})
set(_hdf5_libs hdf5 hdf5_cpp)


find_package(VTK REQUIRED)
vtk_module_config(VTK
  vtkFiltersGeometry #vtkDataSetSurfaceFilter
  vtkFiltersModeling #vtkSphereSource
  vtkIOXML #vtkXMLPolyDataWriter
)
include(${VTK_USE_FILE})

file(GLOB parsers_source_files
    *.cpp
    *.hpp
)

add_library(parsers STATIC ${parsers_source_files})

target_link_libraries(parsers systems translator ${VTK_LIBRARIES} ${_hdf5_libs})

target_include_directories(parsers PUBLIC ${vesicle_SOURCE_DIR}/src)
target_include_directories(parsers PUBLIC ${vesicle_SOURCE_DIR}/analysis/src)
target_include_directories(parsers PUBLIC ${vesicle_SOURCE_DIR}/analysis/HighFive/include/highfive)
