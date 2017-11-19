message(${CMAKE_CURRENT_SOURCE_DIR})

set(IO_source_files
    parameters.cpp
    parameters.hpp
)

add_library(IO STATIC ${IO_source_files})

# target_link_libraries(IO)
# target_include_directories(IO PUBLIC ${PROJECT_SOURCE_DIR}/src)