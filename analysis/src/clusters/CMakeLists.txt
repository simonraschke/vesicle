message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB clusters_source_files
    *.cpp
    *.hpp
)

add_library(clusters STATIC ${clusters_source_files})
target_link_libraries(clusters translator enhance tbb particles)

target_include_directories(clusters PUBLIC ${vesicle_SOURCE_DIR}/src)
target_include_directories(clusters PUBLIC ${vesicle_SOURCE_DIR}/analysis/src)