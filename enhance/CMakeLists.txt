message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB enhance_source_files
    *.cpp
    *.hpp
)

add_library(enhance STATIC ${enhance_source_files})