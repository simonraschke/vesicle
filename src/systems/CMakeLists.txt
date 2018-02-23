message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB systems_source_files
    *.cpp
    *.hpp
)

add_library(systems STATIC ${systems_source_files})

target_link_libraries(systems particles interactions enhance simulations thermostats)
target_include_directories(systems PUBLIC ${PROJECT_SOURCE_DIR}/src)