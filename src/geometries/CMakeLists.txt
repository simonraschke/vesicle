message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB geometries_source_files
    *.cpp
    *.hpp
)

add_library(geometries STATIC ${geometries_source_files})

target_link_libraries(geometries enhance)
target_include_directories(geometries PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(geometries PUBLIC ${PROJECT_SOURCE_DIR}/src)